# ----------------------------------------------------------------------------------------
# simple_recon
# ----------------------------------------------------------------------------------------
# A simple radial reconstruction without SENSE.
#
# Args:
#        rawfile (required)    : The path to the Philips rawfile to be reconstructed
#        output_path (optional): The output path where the results are stored

import argparse
from pathlib import Path

import numpy as np
from image_slide import image_slide
from scipy.io import savemat

import precon as pr

parser = argparse.ArgumentParser(description="normal recon")
parser.add_argument("rawfile", help="path to the raw or lab file")
parser.add_argument(
    "--output-path", default="./", help="path where the output is saved"
)
args = parser.parse_args()

# read parameter
pars = pr.Parameter(Path(args.rawfile))

# define what to read
parameter2read = pr.Parameter2Read(pars.labels)

# enable performance logging (reconstruction times)
pars.performance_logging = True

# dictionary for matlab export
mdic = dict()

# reconstruct every mix and stack seperately
for mix in parameter2read.mix:
    for stack in parameter2read.stack:
        parameter2read.stack = stack
        parameter2read.mix = mix

        # read data
        with open(pars.rawfile, "rb") as raw:
            data, labels = pr.read(raw, parameter2read, pars.labels, pars.coil_info, oversampling_removal=False)

        # sort and zero fill data (create k-space)
        data, labels = pr.sort(data, labels, zeropad=(False, False, False))

        # get the radial trajectory
        data_size = pr.get_data_size(data)
        rot_angle = pars.get_non_cart_angle()
        traj = pr.noncart.get_radial_trajectory(data_size[pr.Enums.X_DIM], data_size[pr.Enums.Y_DIM], rot_angle=rot_angle)

        # Rescale trajectory to match Matlab's normalization (divide by out_res, not no_samples)
        # Matlab: out_res = target_matrix × hardware_ovs (= 240 × 2 = 480 for this dataset)
        target_x = pars.get_recon_resolution(mix=mix, xovs=False)[0]  # base target matrix size (no ovs)
        no_samples = data_size[pr.Enums.X_DIM]  # actual readout length
        out_res_norm = target_x * 2  # = Philips always oversamples radial spokes with factor 2
        traj *= no_samples / out_res_norm  # rescale so that we achieve the desired FOV and resolution after gridding

        # get the radial weights
        weights = pr.noncart.get_radial_weights(traj, data_size[pr.Enums.Y_DIM])

        # set gridder parameters
        output_res = pars.get_recon_resolution(mix=mix, xovs=True, yovs=True, zovs=True)
        osf = 2.0
        kernel_width = 4

        # reshape data to 2D
        data_size = pr.get_data_size(data)
        data = np.reshape(data, (data_size[0]*data_size[1], -1), order='F')

        # grid (returns image directly — no k2i needed)
        data = pr.noncart.nufft(traj, data, weights, output_res, osf, kernel_width)

        # reshape data to nD
        data_size = list(data_size)
        data_size[0] = output_res[0]
        data_size[1] = output_res[1]
        data = np.reshape(data, data_size, order='F')

        # shift data in image space
        yshift = pars.get_shift(enc=1, mix=mix, stack=stack)
        zshift = pars.get_shift(enc=2, mix=mix, stack=stack)
        if yshift:
            data = np.roll(data, yshift, axis=1)
        if zshift:
            data = np.roll(data, zshift, axis=2)

        # partial fourier reconstruction
        kx_range = pars.get_range(enc=0, mix=mix, stack=stack, ovs=False)
        ky_range = pars.get_range(enc=1, mix=mix, stack=stack)
        kz_range = pars.get_range(enc=2, mix=mix, stack=stack)
        if pr.is_partial_fourier(kx_range) or pr.is_partial_fourier(ky_range) or pr.is_partial_fourier(kz_range):
            data = pr.homodyne(data, kx_range, ky_range, kz_range)

        # combine coils with a sum-of squares combination
        data = pr.sos(data, axis=3)

        # perform geometry correction
        r, gys, gxc, gz = pars.get_geo_corr_pars()
        locations = pr.utils.get_unique(labels, "loca")
        MPS_to_XYZ = pars.get_transformation_matrix(loca=locations, mix=mix, target=pr.Enums.XYZ)
        voxel_sizes = pars.get_voxel_sizes(mix=mix)
        data = pr.geo_corr(data, MPS_to_XYZ, r, gys, gxc, gz, voxel_sizes=voxel_sizes)

        mdic[f"data_{mix}_{stack}"] = data
        savemat(Path(args.output_path) / "data.mat", mdic)

        # remove the oversampling
        xovs = pars.get_oversampling(enc=0, mix=mix)
        yovs = pars.get_oversampling(enc=1, mix=mix)
        zovs = pars.get_oversampling(enc=2, mix=mix)
        data = pr.crop(data, axis=(0, 1, 2), factor=(xovs, yovs, zovs), where="symmetric")

        # transform the images into the radiological convention
        transforms = pars.get_in_plane_transformation(mix=mix, stack=stack)
        data = pr.format(data, transforms)

        # make the image square
        res = max(data.shape[0], data.shape[1])
        data = pr.zeropad(data, (res, res), axis=(0, 1))

        # save data in .mat format
        mdic[f"data_{mix}_{stack}"] = data

savemat(Path(args.output_path) / "data.mat", mdic)
