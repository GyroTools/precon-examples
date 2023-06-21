import argparse
from pathlib import Path

import numpy as np
from scipy.io import savemat

import precon as pr
from precon import get_data_size

parser = argparse.ArgumentParser(description='normal recon')
parser.add_argument('rawfile', help='path to the raw or lab file')
parser.add_argument('refscan', help='path to the sense reference scan')
parser.add_argument('--output-path', default='./', help='path where the output is saved')
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

# check if it is a flow scan
segments = parameter2read.extr1
if len(segments) < 2:
    raise RuntimeError('this is not a flow scan')

# dictionary for matlab export
mdic = dict()

# reconstruct every mix and stack seperately
for mix in parameter2read.mix:
    for stack in parameter2read.stack:
        # calculate the sensitivities
        sens = pr.reformat_refscan(qbc, coil, ref_pars, pars, stack=stack, mix=mix, match_target_size=True)

        parameter2read.stack = stack
        parameter2read.mix = mix

        # reconstruct every flow segment separately (to save memory)
        for i in range(0, len(segments)):
            parameter2read.extr1 = segments[i]

            # read data
            with open(args.rawfile, 'rb') as raw:
                data_seg, labels = pr.read(raw, parameter2read, pars.labels, pars.coil_info)

            # sort and zero fill data (create k-space)
            res_before_sense = pars.get_recon_resolution(mix=mix, xovs=False, yovs=True, zovs=True, folded=True)
            data_seg, labels = pr.sort(data_seg, labels, output_size=res_before_sense)

            # ringing filter
            sampled_size = (pars.get_sampled_size(enc=0, stack=stack, ovs=False), pars.get_sampled_size(enc=1, stack=stack),
                            pars.get_sampled_size(enc=2, stack=stack))
            data_seg = pr.hamming_filter(data_seg, (0.25, 0.25, 0.25), axis=(0, 1, 2), sampled_size=sampled_size)

            # FFT
            data_seg = pr.k2i(data_seg, axis=(0, 1, 2))

            # shift data in image space
            yshift = pars.get_shift(enc=1, mix=mix, stack=stack)
            zshift = pars.get_shift(enc=2, mix=mix, stack=stack)
            if yshift:
                data_seg = np.roll(data_seg, yshift, axis=1)
            if zshift:
                data_seg = np.roll(data_seg, zshift, axis=2)

            regularization_factor = pars.goal.get_value(pars.SENSE_REGULARIZATION_FACTOR, at=0, default=2)
            output_size = pars.get_recon_resolution(mix=mix, xovs=False, yovs=True, zovs=True, folded=False)
            data_seg = pr.sense_unfold(data_seg, sens, output_size, regularization_factor=regularization_factor, use_torch=True)

            # partial fourier reconstruction
            kx_range = pars.get_range(enc=0, mix=mix, stack=stack, ovs=False)
            ky_range = pars.get_range(enc=1, mix=mix, stack=stack)
            kz_range = pars.get_range(enc=2, mix=mix, stack=stack)
            if pr.is_partial_fourier(kx_range) or pr.is_partial_fourier(ky_range) or pr.is_partial_fourier(kz_range):
                data_seg = pr.homodyne(data_seg, kx_range, ky_range, kz_range)

            # initialize the final data in the first loop
            if i == 0:
                data_size = list(pr.get_data_size(data_seg))
                data_size[pr.Enums.FLOW_SEGMENT_DIM] = len(segments)
                data = np.zeros(tuple(data_size), dtype=np.csingle, order='F')

            data[:, :, :, :, :, :, :, :, :, [i], ...] = data_seg

        # get the transformation matrices (MPS to XYZ) for every location. it is needed in the geometry correction and
        # the concomitant field correction
        locations = pr.utils.get_unique(labels, 'loca')
        MPS_to_XYZ = pars.get_transformation_matrix(loca=locations, mix=mix, target=pr.Enums.XYZ)
        voxel_sizes = pars.get_voxel_sizes(mix=mix)

        # concommitant field correction (process every location separately)
        concom_factors = pars.get_concom_factors()
        data = pr.concomitant_field_correction(data, MPS_to_XYZ, concom_factors, voxel_sizes, segments)

        # divide the flow segments
        data = pr.divide_flow_segments(data, pars.is_hadamard_encoding())

        # perform geometry correction
        r, gys, gxc, gz = pars.get_geo_corr_pars()
        data = pr.geo_corr(data, MPS_to_XYZ, r, gys, gxc, gz, voxel_sizes=voxel_sizes)

        # remove the oversampling
        yovs = pars.get_oversampling(enc=1, mix=mix)
        zovs = pars.get_oversampling(enc=2, mix=mix)
        data = pr.crop(data, axis=(1, 2), factor=(yovs, zovs), where='symmetric')

        # flow background phase correction
        data = pr.fit_flow_phase(data, order=3)

        # transform the images into the radiological convention
        data = pr.format(data, pars.get_in_plane_transformation(mix=mix, stack=stack))

        # make sure the flow encoding is always along RF-AP-FH axis
        if get_data_size(data, pr.Enums.FLOW_SEGMENT_DIM) <= 3:
            data = pr.format_flow(data, pars.get_coordinate_system(), pars.get_venc(), pars.is_hadamard_encoding())

        # make the image square
        res = max(data.shape[0], data.shape[1])
        data = pr.zeropad(data, (res, res), axis=(0, 1))

        # save data and sensitivities in .mat format
        mdic[f'data_{mix}_{stack}'] = data
        mdic[f'sensitivity_{mix}_{stack}'] = sens.sensitivity
        mdic[f'coil_ref_{mix}_{stack}'] = sens.surfacecoil
        mdic[f'body_ref{mix}_{stack}'] = sens.bodycoil

savemat(Path(args.output_path) / 'data.mat', mdic)