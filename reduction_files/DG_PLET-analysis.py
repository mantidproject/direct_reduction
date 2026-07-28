#================================================================================
#   PLET reduction script v3
#   Uses PolCorr3 routines
#
#   Requires initial LET reduction using 1to1 mapping from DG_reduction.py
#
#                                                                   JRS 18/5/26
#=================================================================================

import sys, os
sys.path.append(os.path.dirname(__file__))
from reduction_utils import PLET_reduce
from mantid.simpleapi import AnalysisDataService as ADS
import time

t = time.time()         # start the clock

#======================================USER INPUT================================
datadir = f'/data/analysis/LET/RBNumber/USER_RB_FOLDER/' # directory where reduced 1to1 .nxs files are stored
sample_runs = range(115331,115339)       # run range
NSF_first = False                        # this flag should be true if you measured with the flipper off first (check JournalViewer)
pressure  = 0.9                          # the pressure of the 3He gas
he_mode   = 'fit'                        # options: 'fit' or 'direct' -
                                         # use 'fit' for set_helium_parameters
eis       = [8.61, 3.60, 1.97, 1.24]     # incident energies as a list (cannot use "auto")
PF        = [0.90, 0.91, 0.92, 0.93]     # Polarizer*Flipper efficiency (one for each Ei)
rings_map = 'RINGS_MAP_XML'              # Set to None for 1to1 output - must be a .map file
input_format = '.nxspe'                  # format of input files
                                         # if '.nxspe' then LET_rings_153.map will be used for rings grouping
output_format = '.nxspe'                 # format of output files - either '.nxspe' or '.nxs'
output_label  = 'quartz'                 # file stem label for output files
calib = "3HeCal_115331-115338.txt"       # cell calibration file (made using get_helium_parameters())
cycle = 'CYCLE_ID'                       # cycle number
#=================================END OF USER INPUT==============================

# Add necessary folders to Mantid path if they have not been added
cycle_shortform = cycle[2:] if cycle.startswith('20') else cycle
config.appendDataSearchDir(f'/archive/NDXLET/Instrument/data/cycle_{cycle_shortform}/')
config.appendDataSearchDir(f'/data/instrument/{inst}/CYCLE20{cycle_shortform.replace("_","")}/USER_RB_FOLDER/')
config.appendDataSearchDir(f'/usr/local/mprogs/InstrumentFiles/{inst.swapcase()}/')

# The next line clears *all* workspaces ("temporary" fix for memory leak), which will
# cause issues for MSlice: you should run MSlice in a separate Mantid process to this script.
ADS.clear()
    
data = PLET_reduce(sample_runs, eis, PF,
                   he_pressure = pressure,
                   he_mode     = he_mode,
                   NSF_first   = NSF_first,
                   label       = output_label,
                   rings_map   = rings_map,
                   datadir     = datadir,
                   file_format = input_format)
 
#data.get_PF_from_monitor(PHe=0.58)      # Only for instrument scientist
#data.get_PF_from_quartz()               # Only for instrument scientist
data.get_helium_parameters(save=True)    # This saves a "3HeCal_<first_run>-<last_run>.txt" file
#data.set_helium_parameters(PHe0=0.58, T1=60, T0run=115331)  # Only for instrument scientist
#data.set_helium_parameters(cal=calib)   # This loads a previously calculated calibration file
data.correct_data()
data.components()
data.output(type=output_format)

#=====================================================================================

print(f'\nPLET-analysis: Reduction complete in {time.time() - t:.1f} seconds\n')
