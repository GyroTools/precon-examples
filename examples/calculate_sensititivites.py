# ----------------------------------------------------------------------------------------
# calculate_sensitivities
# ----------------------------------------------------------------------------------------
# Calculates the sensitivity maps from a Philips reference scan
#
# Args:
#        refscan (required)          : The path to the Philips SENSE reference scan
#        target_scan (required)      : The path to the SENSE scan (target scan)
#        output_path (optional)      : The output path where the results are stored
#        match-target-size (optional): When given the sensitivity maps have the same size as the target images
#        fov (optional)              : A user defined field-of-view for the sensitivity maps
#        output_size (optional)      : A user defined output size for the sensitivity maps

import argparse
from pathlib import Path

from scipy.io import savemat

from precon import calculate_sensitivities

parser = argparse.ArgumentParser(description='normal recon')
parser.add_argument('refscan', help='path to the raw or lab file of the sense refscan')
parser.add_argument('target_scan', help='path to the raw or lab file of the target scan')
parser.add_argument('--output-path', default='./', help='path where the output is saved')
parser.add_argument('--match-target-size', action='store_true', help='if given then the size of the sensitivities matches the one of the target data')
parser.add_argument('--fov', nargs="+", type=float, default=None, help='the FOV of the sensitivity maps')
parser.add_argument('--output_size', nargs="+", type=int, default=None, help='the putput size of the sensitivities')
args = parser.parse_args()

s = calculate_sensitivities(Path(args.refscan), Path(args.target_scan), stack=0, mix=0, match_target_size=args.match_target_size, fov=args.fov, output_size=args.output_size)

mdic = {'qbc': s.bodycoil, 'coil': s.surfacecoil, 'sensitivity': s.sensitivity, 'psi': s.psi}
savemat(Path(args.output_path) / 'sense.mat', mdic)