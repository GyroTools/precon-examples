# ----------------------------------------------------------------------------------------
# cardiac_retro_recon
# ----------------------------------------------------------------------------------------
# Performs a reconstruction of a retrospectively triggered cardiac acquisition
#
# Args:
#        rawfile (required)    : The path to the Philips rawfile to be reconstructed
#        output_path (optional): The output path where the results are stored
#        nr_phases   (optional): The number of cardiac phases to be reconstructed
#
# The reconstruction performed in this file consists of the following steps:
#
#   1. Read the parameters from the sinfile
#   2. Create a Parameter2Read class from the labels which defines what data to read
#   3. Loop over all mixes and stacks
#   4. Read the data from the current mix and stack (the basic corrections as well as the oversampling removal in readout direction is performed in the reader)
#   5. Cardiac binning. A heart phase number is assigned to every acquired k-space profile.
#   6. Sort and zero-fill the data according to the labels (create k-space)
#   7. Fill the holes k-space due to an irregular heart rate, or if the reconstructed heart phases are larger than the acquired ones. The holes are filled by a nearest neighbour interpolation over time.
#   8. Perform fourier transformation
#   9. Shift the images such that they are aligned correctly
#  10. Perform a partial fourier (homodyne) reconstruction when halfscan or partial echo was enabled
#  11. Combine the coils with a sum-of-squares combination
#  12. Perform the geometry correction
#  13. Remove the oversampling along the phase encoding directions
#  14. Transform the images into the radiological convention
#  15. Make the images square

import argparse
from pathlib import Path

import numpy as np
from scipy.io import savemat

import precon as pr


parser = argparse.ArgumentParser(description='normal recon')
parser.add_argument('rawfile', help='path to the raw or lab file')
parser.add_argument('--output-path', default='./', help='path where the output is saved')
parser.add_argument('--nr_phases', type=int, default=None, help='the number of cardiac phases to be reconstructed')
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
        with open(args.rawfile, 'rb') as raw:
            data, labels = pr.read(raw, parameter2read, pars.labels, pars.coil_info)

        # retrospective cardiac binning
        nr_phases = args.nr_phases if args.nr_phases else pars.get_nr_phases()
        labels = pr.retro_binning(labels, nr_phases)

        # sort and zero fill data (create k-space)
        cur_recon_resolution = pars.get_recon_resolution(mix=mix, xovs=False, yovs=True, zovs=True)
        data, labels = pr.sort(data, labels, output_size=cur_recon_resolution)

        # fill the holes in k-space due to retrospective binning
        data = pr.retro_fill_holes(data)

        mdic[f'data'] = data
        savemat(Path(args.output_path) / 'kspace.mat', mdic)

        # FFT
        data = pr.k2i(data, axis=(0, 1, 2))

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
        locations = pr.utils.get_unique(labels, 'loca')
        MPS_to_XYZ = pars.get_transformation_matrix(loca=locations, mix=mix, target=pr.Enums.XYZ)
        voxel_sizes = pars.get_voxel_sizes(mix=mix)
        data = pr.geo_corr(data, MPS_to_XYZ, r, gys, gxc, gz, voxel_sizes=voxel_sizes)

        # remove the oversampling
        yovs = pars.get_oversampling(enc=1, mix=mix)
        zovs = pars.get_oversampling(enc=2, mix=mix)
        data = pr.crop(data, axis=(1, 2), factor=(yovs, zovs), where='symmetric')

        # transform the images into the radiological convention
        data = pr.format(data, pars.get_in_plane_transformation(mix=mix, stack=stack))

        # make the image square
        res = max(data.shape[0], data.shape[1])
        data = pr.zeropad(data, (res, res), axis=(0, 1))

        # save data in .mat format
        mdic[f'data_{mix}_{stack}'] = data

savemat(Path(args.output_path) / 'data.mat', mdic)
