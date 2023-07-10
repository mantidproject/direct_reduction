from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import numpy as np
import re, os, h5py

#========================================================
# Functions to copy instrument info needed by HORACE
# Resolution convolution if it exists in the raw file
def get_raw_file_from_ws(ws):
    for alg in [h for h in ws.getHistory().getAlgorithmHistories() if 'Load' in h.name()]:
        for prp in [a for a in alg.getProperties() if 'Filename' in a.name()]:
            if re.search('[0-9]*.nxs', prp.value()) is not None:
                return prp.value()
    raise RuntimeError('Could not find raw NeXus file in workspace history')

def copy_inst_info(outfile, in_ws):
    try:
        raw_file_name = get_raw_file_from_ws(mtd[in_ws])
    except RuntimeError:
        return
    print(raw_file_name)
    if not os.path.exists(outfile):
        outfile = os.path.join(mantid.simpleapi.config['defaultsave.directory'], os.path.basename(outfile))
    print(outfile)
    with h5py.File(raw_file_name, 'r') as raw:
        exclude = ['dae', 'detector_1', 'name']
        to_copy = [k for k in raw['/raw_data_1/instrument'] if not any([x in k for x in exclude])]
        if 'aperture' not in to_copy and 'mono_chopper' not in to_copy:
            return
        with h5py.File(outfile, 'r+') as spe:
            print(spe.keys())
            spe_root = list(spe.keys())[0]
            en0 = spe[f'{spe_root}/instrument/fermi/energy'][()]
            if 'fermi' in to_copy:
                del spe[f'{spe_root}/instrument/fermi']
            for grp in to_copy:
                print(grp)
                src = raw[f'/raw_data_1/instrument/{grp}']
                h5py.Group.copy(src, src, spe[f'{spe_root}/instrument/'])
            if 'fermi' in to_copy:
                spe[f'{spe_root}/instrument/fermi/energy'][()] = en0
            detroot = f'{spe_root}/instrument/detector_elements_1'
            spe.create_group(detroot)
            for df0, df1 in zip(['SPEC', 'UDET', 'DELT', 'LEN2', 'CODE', 'TTHE', 'UT01'], \
                ['spectrum_number', 'detector_number', 'delt', 'distance', 'detector_code', 'polar_angle', 'azimuthal_angle']):
                src = raw[f'/raw_data_1/isis_vms_compat/{df0}']
                h5py.Group.copy(src, src, spe[detroot], df1)
            for nn in range(raw['/raw_data_1/isis_vms_compat/NUSE'][0]):
                src = raw[f'/raw_data_1/isis_vms_compat/UT{nn+1:02d}']
                h5py.Group.copy(src, src, spe[detroot], f'user_table{nn+1:02d}')
#========================================================

