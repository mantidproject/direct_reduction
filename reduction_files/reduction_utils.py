from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import numpy as np
import re, os, h5py

#========================================================
# General utility functions

def rename_existing_ws(ws_name):
    # Gets an existing workspace and ensures we have its monitors for reduction
    ws = CloneWorkspace(ws_name, OutputWorkspace='ws')
    try:
        mon_ws = ws.getMonitorWorkspace()
    except RuntimeError:
        specInfo = ws.spectrumInfo()
        try:
            mon_list = [i for i in range(ws.getNumberHistograms()) if specInfo.isMonitor(i)]
        except RuntimeError:
            mon_list = []
        if len(mon_list) > 0:
            ExtractMonitors(ws_name, DetectorWorkspace='ws', MonitorWorkspace='ws_monitors')
        else:
            _get_mon_from_history(ws_name)
    else:
        CloneWorkspace(mon_ws, OutputWorkspace='ws_monitors')

def _get_mon_from_history(ws_name):
    # Tries to look in the history of a workspace for its original raw file
    # loads the monitors from that raw file.
    orig_file = None
    for hist in mtd[ws_name].getHistory().getAlgorithmHistories():
        if hist.name().startswith('Load') and 'Filename' in [pp.name() for pp in hist.getProperties()]:
            orig_file = hist.getPropertyValue('Filename')
            break
    if orig_file is None:
        raise RuntimeError(f'Cannot find original file from workspace {ws_name} to load logs from')
    try:
        Load(orig_file, SpectrumMax=10, LoadMonitors=True, OutputWorkspace='tmp_mons')
    except TypeError:
        Load(orig_file, SpectrumMax=10, LoadMonitors='Separate', OutputWorkspace='tmp_mons')
    ws_mon_name = f'{ws_name}_monitors'
    RenameWorkspace('tmp_mons_monitors', ws_mon_name)
    DeleteWorkspace('tmp_mons')
    mtd[ws_name].setMonitorWorkspace(mtd[ws_mon_name])
    CloneWorkspace(ws_mon_name, OutputWorkspace='ws_monitors')

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
# MARI specific functions

def remove_extra_spectra_if_mari(wsname='ws'):
    ws = mtd[wsname]
    if ws.getNumberHistograms() > 918:
        ws = RemoveSpectra(ws, [0])
    try:
        ws_mon = ws.getMonitorWorkspace()
    except RuntimeError:
        ws_mon = mtd[f'{wsname}_monitors']
    if ws_mon.getNumberHistograms() > 3:
        ws_monitors = RemoveSpectra(ws_mon, [3])

def shift_frame_for_mari_lowE(origEi, wsname='ws_norm', wsmon='ws_monitors'):
    ws_norm, ws_monitors = (mtd[wsname], mtd[wsmon])
    if origEi < 4.01:
        # If Ei < 4, mon 3 is in 2nd frame so need to shift it by 20ms
        ws_monitors = ScaleX(wsmon, 20000, Operation='Add', IndexMin=2, IndexMax=2, OutputWorkspace=wsmon)
        if origEi < 3.1:
            # Additionally if Ei<3.1, data will also be in 2nd frame, shift all ToF by 20ms
            ws_norm = ScaleX(wsname, 20000, Operation='Add', IndexMin=0, IndexMax=ws_norm.getNumberHistograms()-1, OutputWorkspace=wsname)
    return ws_norm, ws_monitors

#========================================================
