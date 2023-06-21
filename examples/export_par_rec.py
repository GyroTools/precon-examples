from pathlib import Path

import numpy as np
from scipy.io import savemat

import precon as pr
from precon import get_data_size

folder = Path(r'C:\Users\martin\Desktop\temp\work')
data = np.load(str(folder / 'data.npy'), allow_pickle=True)
labels = np.load(str(folder / 'labels.npy'), allow_pickle=True)
venc = np.load(str(folder / 'venc.npy'), allow_pickle=True)

pars = pr.Parameter(Path(r'C:\Users\martin\Desktop\temp\csflow\phantom\cm_11072022_1648068_11_1_wip_5d-qfl-3.5xV4.raw'))

# --------------
# Par/Rec Export
print('Par/Rec Export')
types = (pr.Enums.REC_IMAGE_TYPE_M, pr.Enums.REC_IMAGE_TYPE_P)
basename = pars.rawfile.stem

scaling = pr.Recfile.get_scaling(data, types=types)
for s in range(0, get_data_size(data, pr.Enums.FLOW_SEGMENT_DIM)):
    rec = pr.Recfile(data[:, :, :, :, :, :, :, :, :, s, :, :], types=types, scaling=scaling)
    savemat(folder / 'data_rec.mat', {'data_rec': rec})
    par = pr.Parfile(pars, data[:, :, :, :, :, :, :, :, :, [s], :, :], labels[:, :, :, :, :, :, :, :, :, [s], :, :], types=types, scaling=scaling)
    par.General.PhaseEncodingVelocity = venc[s,:]
    filename_par = folder / f'{basename}_{s}.par'
    filename_rec = folder / f'{basename}_{s}.rec'
    par.write(filename_par)
    rec.write(filename_rec)