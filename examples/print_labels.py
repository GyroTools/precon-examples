import argparse
from pathlib import Path

from prettytable import PrettyTable

import precon as pr

parser = argparse.ArgumentParser(description="normal recon")
parser.add_argument("rawfile", help="path to the raw or lab file")

args = parser.parse_args()

# read parameter
pars = pr.Parameter(Path(args.rawfile))

labels = pars.labels

# convert labels to list of dicts
labels_dict = [{field[0]: getattr(label, field[0]) for field in label._fields_} for label in labels]

table = PrettyTable()

headers = labels_dict[0].keys()
table.field_names = headers

# Add rows to the table
for label in labels_dict:
    table.add_row(label.values())

# Print the table
print(table)