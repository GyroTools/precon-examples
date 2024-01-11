# ----------------------------------------------------------------------------------------
# spectro_recon
# ----------------------------------------------------------------------------------------
# A simple reconstruction for single-voxel spectroscopy data
#
# Args:
#        rawfile (required)    : The path to the Philips rawfile to be reconstructed
#        output_path (optional): The output path where the results are stored

import argparse
from pathlib import Path

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

# reconstruct the first mix and stack
parameter2read.stack = 0
parameter2read.mix = 0

# read data
with open(pars.rawfile, "rb") as raw:
    data, labels = pr.read(raw, parameter2read, pars.labels, pars.coil_info, oversampling_removal=False)

data, labels = pr.sort(data, labels, immediate_averaging=True, zeropad=(False, False, False))

# determine the oversampling factor and remove oversampling
cur_recon_resolution = pars.get_recon_resolution()
xovs = pr.get_data_size(data)[pr.Enums.X_DIM] / cur_recon_resolution[pr.Enums.X_DIM]
data = pr.spectro_downsample(data, xovs)

# combine coils with a svd combination
data = pr.spectro_combine_coils(data)

# FFT
data = pr.k2i(data, axis=(0, 1, 2))

# save data in .mat format
mdic[f"data"] = data

savemat(Path(args.output_path) / "data.mat", mdic)
