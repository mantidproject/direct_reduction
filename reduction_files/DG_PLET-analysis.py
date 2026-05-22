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
datadir    = f'/home/rs5329/Documents/Team_LET/PLET_quartz/'   #directory where reduced 1to1 .nxs files are stored
NSF_first = False                        # this flag should be true if you measured with the flipper off first (check JournalViewer)
pressure  = 0.9                          # the pressure of the 3He gas
he_mode   = 'fit'                        # options: 'fit' or 'direct' -
                                         # use 'fit' for set_helium_parameters
eis       = [8.61, 3.60, 1.97, 1.24]     # incident energies
PF        = [0.90, 0.91, 0.92, 0.93]     # Polarizer*Flipper efficiency   
rings_map = 'LET_rings_261.xml'          # Set to None for 1to1 output - must be a .map file
format    = '.nxspe'                     # format of input files
                                         # if '.nxspe' then LET_rings_153.map will be used for rings grouping
output_format = '.nxspe'                 # format of output files - either '.nxspe' or '.nxs'
output_label  = 'quartz'                # file stem label for output files
cell_calib = "3HeCal_115331-115412.txt"  # cell calibration file (made using get_helium_parameters())
#PHe0=0.6; T1=60.; T0run=115331          # cell parameters input directly
sample_runs = range(115331,115339)       # run range
#=================================END OF USER INPUT===================================

ADS.clear() # This clear all the workspaces ("temporary" fix for memory leak)
    
data = PLET_reduce(sample_runs, eis, PF,
                   he_pressure = pressure,
                   he_mode     = he_mode,
                   NSF_first   = NSF_first,
                   label       = output_label,
                   rings_map   = rings_map,
                   datadir     = datadir,
                   file_format = format)
                           
#data.get_PF_from_monitor(PHe0=0.58)
#data.get_PF_from_quartz()
#data.get_helium_parameters()
#data.set_helium_parameters(PHe0=0.58,T1=60,T0run=115331)
data.set_helium_parameters(cal=cell_calib)
data.correct_data()
data.components()
data.output(type=output_format)

#=====================================================================================

print(f'\nPLET-analysis: Reduction complete in {time.time() - t:.1f} seconds\n')
