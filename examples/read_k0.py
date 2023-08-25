# ----------------------------------------------------------------------------------------
# read_k0
# ----------------------------------------------------------------------------------------
# Reads all central k-space profiles (ky = kz = 0)
#
# Args:
#        rawfile (required)    : The path to the Philips rawfile to be reconstructed

import argparse
from pathlib import Path

import precon as pr

parser = argparse.ArgumentParser(description='normal recon')
parser.add_argument('rawfile', help='path to the raw or lab file')
args = parser.parse_args()

# read parameter
pars = pr.Parameter(Path(args.rawfile))

labels = pars.labels.copy()

# set the labels of all rejected profiles to the normal type
for label in labels:
    if label.typ == pr.Enums.REJECTED_DATA:
        label.typ = pr.Enums.NORMAL_DATA

parameter2read = pr.Parameter2Read(labels)
parameter2read.typ = pr.Enums.NORMAL_DATA
parameter2read.ky = 0
parameter2read.kz = 0

# read central k-space lines
with open(pars.rawfile, 'rb') as raw:
    data, labels = pr.read(raw, parameter2read, labels, pars.coil_info)
