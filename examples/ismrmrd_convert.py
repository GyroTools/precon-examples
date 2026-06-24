# ----------------------------------------------------------------------------------------
# ismrmrd_convert
# ----------------------------------------------------------------------------------------
# Convert a Philips raw MRI dataset to the ISMRMRD HDF5 format.
#
# Args:
#        rawfile (required)     : The path to the Philips rawfile to be converted
#        output_file (optional) : The output .h5 file path (default: <rawfile>.h5)
#        --no-anonymize         : Include patient name, ID and birthdate in the header
#        --no-overwrite         : Raise an error if the output file already exists
#        --block-size-mb        : Read chunk size in MB (default: 100)

import argparse
from pathlib import Path

from precon.convert.ismrmrd_converter import ismrmrd_converter

parser = argparse.ArgumentParser(description="Convert Philips raw data to ISMRMRD HDF5")
parser.add_argument("rawfile", help="path to the Philips .raw file")
parser.add_argument(
    "output_file",
    nargs="?",
    default=None,
    help="output .h5 file path (default: <rawfile>.h5)",
)
parser.add_argument(
    "--no-anonymize",
    action="store_true",
    help="include patient name, ID and birthdate in the header",
)
parser.add_argument(
    "--no-overwrite",
    action="store_true",
    help="raise an error if the output file already exists",
)
parser.add_argument(
    "--block-size-mb",
    type=int,
    default=100,
    help="read chunk size in MB (default: 100)",
)
args = parser.parse_args()

rawfile = Path(args.rawfile)
output_file = Path(args.output_file) if args.output_file else rawfile.with_suffix(".h5")

conv = ismrmrd_converter(
    rawfile=rawfile,
    block_size_mb=args.block_size_mb,
    anonymize=not args.no_anonymize,
    overwrite=not args.no_overwrite,
)
conv.convert(output_file)
print(f"Written: {output_file}")
