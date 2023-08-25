# ----------------------------------------------------------------------------------------
# epi_recon
# ----------------------------------------------------------------------------------------
# An EPI reconstruction with optional SENSE reconstruction
#
# Args:
#        rawfile (required)    : The path to the Philips rawfile to be reconstructed
#        refscan (optional)    : The path to the Philips SENSE reference scan
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
#   7. Read the epi correction data
#   8. Grid the samples on a regular grid using the non-uniform-sampling coordinates (needed because data was sampled on the gradient ramp)
#   9. Sort and zero-fill the data according to the labels (create k-space)
#  10. Sort the epi correction data
#  11. Perform fourier transformation along readout direction
#  12. Shift the data along the readout direction
#  13. Perform the epi correction
#  14. Perform fourier transformation along phase encoding directions
#  15. Shift the data along phase encoding directions
#  16. Perform a SENSE reconstruction (unfolding)
#  17. Perform a partial fourier (homodyne) reconstruction when halfscan or partial echo was enabled
#  18. Perform the geometry correction
#  19. Remove the oversampling along the phase encoding directions
#  20. Transform the images into the radiological convention
#  21. Make the images square

import argparse
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.io import savemat

import precon as pr

parser = argparse.ArgumentParser(description='normal recon')
parser.add_argument('rawfile', help='path to the raw or lab file')
parser.add_argument('--refscan', default=None, help='path to the sense reference scan')
parser.add_argument('--output-path', default='./', help='path where the output is saved')
args = parser.parse_args()

# read parameter
pars = pr.Parameter(Path(args.rawfile))

# define what to read
parameter2read = pr.Parameter2Read(pars.labels)

# enable performance logging (reconstruction times)
pars.performance_logging = True

if args.refscan:
    # reconstruct refscan
    ref_pars = pr.Parameter(Path(args.refscan))
    qbc, coil = pr.reconstruct_refscan(ref_pars)

# check if it is an epi scan
if not pars.is_epi():
    raise RuntimeError('this is not an epi scan')

# dictionary for matlab export
mdic = dict()

sens = None
sense_factors = None

# reconstruct every mix and stack seperately
for mix in parameter2read.mix:
    for stack in parameter2read.stack:
        if args.refscan:
            # calculate the sensitivities
            sens = pr.reformat_refscan(qbc, coil, ref_pars, pars, stack=stack, mix=mix, match_target_size=True)
            sense_factors = pars.get_value(pars.SENSE_FACTORS, default=[1, 1, 1])

        parameter2read.stack = stack
        parameter2read.mix = mix

        # read data
        with open(args.rawfile, 'rb') as raw:
            data, labels = pr.read(raw, parameter2read, pars.labels, pars.coil_info, oversampling_removal=False)

            # read epi correction data
            parameter2read.typ = pr.Label.TYPE_ECHO_PHASE
            epi_corr_data, epi_corr_labels = pr.read(raw, parameter2read, pars.labels, pars.coil_info, oversampling_removal=False)

        # grid the data from the nus encoding numbers to a regular grid
        nus_enc_nrs = pars.get_nus_enc_nrs()
        kx_range = pars.get_range(mix=mix, stack=stack)
        f = interp1d(nus_enc_nrs, data, axis=0, bounds_error=False, fill_value=0)
        data = f(np.arange(kx_range[0], kx_range[1]+1))
        f = interp1d(nus_enc_nrs, epi_corr_data, axis=0, bounds_error=False, fill_value=0)
        epi_corr_data = f(np.arange(kx_range[0], kx_range[1] + 1))


        # sort and zero fill data (create k-space)
        cur_recon_resolution = pars.get_recon_resolution(mix=mix, xovs=True, yovs=True, zovs=True)
        data, labels = pr.sort(data, labels, output_size=cur_recon_resolution)

        # sort the epi correction data (since ky is always 0 set the grad label as ky)
        epi_corr_data, epi_corr_labels = pr.sort(epi_corr_data, epi_corr_labels, output_size=[cur_recon_resolution[0]], zeropad=(True, False, False), immediate_averaging=False, ky='grad')

        # FFT along readout direction
        data = pr.k2i(data, axis=0)
        epi_corr_data = pr.k2i(epi_corr_data, axis=0)

        # shift data in image space
        xshift = pars.get_shift(enc=0, mix=mix, stack=stack)
        if xshift:
            data = np.roll(data, xshift, axis=0)
            epi_corr_data = np.roll(epi_corr_data, xshift, axis=0)

        # epi correction
        slopes, offsets = pr.get_epi_corr_data(epi_corr_data, epi_corr_labels)
        data = pr.epi_corr(data, labels, slopes, offsets)

        # FFT along phase encoding direction
        data = pr.k2i(data, axis=(1, 2))

        # shift data in image space
        yshift = pars.get_shift(enc=1, mix=mix, stack=stack)
        zshift = pars.get_shift(enc=2, mix=mix, stack=stack)
        if yshift:
            data = np.roll(data, yshift, axis=1)
        if zshift:
            data = np.roll(data, zshift, axis=2)

        # remove the oversampling along readout direction
        xovs = pars.get_oversampling(enc=0, mix=mix)
        data = pr.crop(data, axis=0, factor=xovs, where='symmetric')

        # SENSE unfold
        if args.refscan:
            regularization_factor = pars.get_value(pars.SENSE_REGULARIZATION_FACTOR, at=0, default=2)
            data = pr.sense_unfold(data, sens, sense_factors, regularization_factor=regularization_factor)

        # partial fourier reconstruction
        kx_range = pars.get_range(enc=0, mix=mix, stack=stack, ovs=False)
        ky_range = pars.get_range(enc=1, mix=mix, stack=stack)
        kz_range = pars.get_range(enc=2, mix=mix, stack=stack)
        if pr.is_partial_fourier(kx_range) or pr.is_partial_fourier(ky_range) or pr.is_partial_fourier(kz_range):
            data = pr.homodyne(data, kx_range, ky_range, kz_range)

        # combine coils with a sum-of squares combination
        if not args.refscan:
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
        data = pr.zeropad(data, (res, res), axis=(0,1))

        # save data in .mat format
        mdic[f'data_{mix}_{stack}'] = data

savemat(Path(args.output_path) / 'data.mat', mdic)