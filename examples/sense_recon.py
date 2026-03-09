# ----------------------------------------------------------------------------------------
# sense_recon
# ----------------------------------------------------------------------------------------
# A cartesian SENSE (or CLEAR) reconstruction which calculates the sensitivity maps from the
# rawfile of the SENSE reference scan
#
# Args:
#        rawfile (required)    : The path to the Philips rawfile to be reconstructed
#        refscan (required)    : The path to the Philips SENSE reference scan
#        output_path (optional): The output path where the results are stored
#
# The reconstruction performed in this file consists of the following steps:
#
#   1. Read the parameters from the rawfile
#   2. Reconstruct the SENSE reference scan
#   3. Create a Parameter2Read class from the labels which defines what data to read
#   4. Loop over all mixes and stacks
#   5. Reformat the SENSE reference scan into the geometry of the target scan
#   6. Read the data from the current mix and stack (the basic corrections as well as the oversampling removal in readout direction is performed in the reader)
#   7. Sort and zero-fill the data according to the labels (create k-space)
#   8. Apply a ringing filter
#   9. Perform fourier transformation
#  10. Shift the images such that they are aligned correctly
#  11. Perform a SENSE reconstruction (unfolding)
#  12. Perform a partial fourier (homodyne) reconstruction when halfscan or partial echo was enabled
#  13. Perform the geometry correction
#  14. Remove the oversampling along the phase encoding directions
#  15. Transform the images into the radiological convention
#  16. Make the images square

import argparse
from pathlib import Path

import numpy as np
from scipy.io import savemat

import precon as pr

parser = argparse.ArgumentParser(description="normal recon")
parser.add_argument("rawfile", help="path to the raw or lab file")
parser.add_argument("refscan", help="path to the sense reference scan")
parser.add_argument(
    "--output-path", default="./", help="path where the output is saved"
)
args = parser.parse_args()

# read parameter
pars = pr.Parameter(Path(args.rawfile))

# enable performance logging (reconstruction times)
pars.performance_logging = True

# reconstruct refscan
ref_pars = pr.Parameter(Path(args.refscan))
qbc, coil = pr.reconstruct_refscan(ref_pars)

# define what to read
parameter2read = pr.Parameter2Read(pars.labels)

# dictionary for matlab export
mdic = dict()

# reconstruct every mix and stack seperately
for mix in parameter2read.mix:
    for stack in parameter2read.stack:
        parameter2read.stack = stack
        parameter2read.mix = mix

        # calculate the sensitivities
        sens = pr.reformat_refscan(qbc, coil, ref_pars, pars, stack=stack, mix=mix, match_target_size=True)

        # read data
        with open(pars.rawfile, "rb") as raw:
            data, labels = pr.read(raw, parameter2read, pars.labels, pars.coil_info)

        # sort and zero fill data (create k-space)
        res_before_sense = pars.get_recon_resolution(mix=mix, xovs=False, yovs=True, zovs=True, folded=True)
        data, labels = pr.sort(data, labels, output_size=res_before_sense)

        # ringing filter
        sampled_size = (pars.get_sampled_size(enc=0, stack=stack, ovs=False), pars.get_sampled_size(enc=1, stack=stack),
                        pars.get_sampled_size(enc=2, stack=stack))
        data = pr.hamming_filter(data, (0.25, 0.25, 0.25), axis=(0, 1, 2), sampled_size=sampled_size)

        # FFT
        data = pr.k2i(data, axis=(0, 1, 2))

        # shift data in image space
        yshift = pars.get_shift(enc=1, mix=mix, stack=stack)
        zshift = pars.get_shift(enc=2, mix=mix, stack=stack)
        if yshift:
            data = np.roll(data, yshift, axis=1)
        if zshift:
            data = np.roll(data, zshift, axis=2)

        # SENSE unfolding
        regularization_factor = pars.get_value(pars.SENSE_REGULARIZATION_FACTOR, at=0, default=2)
        output_size = pars.get_recon_resolution(mix=mix, xovs=False, yovs=True, zovs=True, folded=False)
        data = pr.sense_unfold(data, sens, output_size, regularization_factor=regularization_factor)

        # partial fourier reconstruction
        kx_range = pars.get_range(enc=0, mix=mix, stack=stack, unfolded=True, ovs=False)
        ky_range = pars.get_range(enc=1, mix=mix, stack=stack, unfolded=True)
        kz_range = pars.get_range(enc=2, mix=mix, stack=stack, unfolded=True)
        if pr.is_partial_fourier(kx_range) or pr.is_partial_fourier(ky_range) or pr.is_partial_fourier(kz_range):
            data = pr.homodyne(data, kx_range, ky_range, kz_range)

        # perform geometry correction
        r, gys, gxc, gz = pars.get_geo_corr_pars()
        locations = pr.utils.get_unique(labels, "loca")
        MPS_to_XYZ = pars.get_transformation_matrix(loca=locations, mix=mix, target=pr.Enums.XYZ)
        voxel_sizes = pars.get_voxel_sizes(mix=mix)
        data = pr.geo_corr(data, MPS_to_XYZ, r, gys, gxc, gz, voxel_sizes=voxel_sizes)

        # remove the oversampling
        yovs = pars.get_oversampling(enc=1, mix=mix)
        zovs = pars.get_oversampling(enc=2, mix=mix)
        data = pr.crop(data, axis=(1, 2), factor=(yovs, zovs), where="symmetric")

        # transform the images into the radiological convention
        data = pr.format(data, pars.get_in_plane_transformation(mix=mix, stack=stack))

        # make the image square
        res = max(data.shape[0], data.shape[1])
        data = pr.zeropad(data, (res, res), axis=(0, 1))

        # save data and sensitivities in .mat format
        mdic[f"data_{mix}_{stack}"] = data
        mdic[f"sensitivity_{mix}_{stack}"] = sens.sensitivity
        mdic[f"coil_ref_{mix}_{stack}"] = sens.surfacecoil
        mdic[f"body_ref{mix}_{stack}"] = sens.bodycoil

savemat(Path(args.output_path) / "data.mat", mdic)
