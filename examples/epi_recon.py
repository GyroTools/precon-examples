import argparse
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.io import savemat

import precon as pr
from precon.parameter.models.label import Label

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
            sense_factors = pars.goal.get_value(pars.SENSE_FACTORS, default=[1, 1, 1])

        parameter2read.stack = stack
        parameter2read.mix = mix

        # read data
        with open(args.rawfile, 'rb') as raw:
            data, labels = pr.read(raw, parameter2read, pars.labels, pars.coil_info, oversampling_removal=False)

            # read epi correction data
            parameter2read.typ = Label.TYPE_ECHO_PHASE
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
            regularization_factor = pars.goal.get_value(pars.SENSE_REGULARIZATION_FACTOR, at=0, default=2)
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