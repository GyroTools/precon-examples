import argparse
from pathlib import Path

from scipy.io import savemat

from precon import calculate_sensitivities

parser = argparse.ArgumentParser(description='normal recon')
parser.add_argument('refscan', help='path to the raw or lab file of the sense refscan')
parser.add_argument('target_scan', help='path to the raw or lab file of the target scan')
parser.add_argument('--output-path', default='./', help='path where the output is saved')
parser.add_argument('--match-target-size', action='store_true', help='if given then the size of the sensitivities matches the one of the target data')
args = parser.parse_args()

s = calculate_sensitivities(Path(args.refscan), Path(args.target_scan), stack=0, mix=0, match_target_size=args.match_target_size)

mdic = {'qbc': s.bodycoil, 'coil': s.surfacecoil, 'sensitivity': s.sensitivity, 'psi': s.psi}
savemat(Path(args.output_path) / 'sense.mat', mdic)