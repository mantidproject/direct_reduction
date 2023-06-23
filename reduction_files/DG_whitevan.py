# =======================================================
#
#      WHITE VAN REDUCTION SCRIPT
#      JRS 26/4/23
#
#      Reads white van run, subtracts empty, and
#      writes out white vanadium integrals to an
#      ascii file for use with DG_reduction.py
#
#      No mask is employed at this stage
#
#========================================================

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import matplotlib.pyplot as plt
import numpy as np
import time

t = time.time()         #start the clock

#=======================User Inputs======================
whitevan       = 78715                  #white vanadium run
whitevan_bg    = None                   #background for the white vanadium
whitevan_trans = 1                      #transmission factor
#========================================================

#==================Local Contact Inputs==================
config['default.instrument']='LET'   #instrument
cycle = '22_1'                          #cycle
wv_lrange = [1,5]                     #wavelength integration limits for output of vanadium integrals
wv_detrange = [30000,60000]             #spectrum index range for average intensity calculation
idebug = False                          #keep itermediate workspaces for debugging
#========================================================

inst = config['default.instrument']
datadir = '/archive/NDX'+inst+'/instrument/data/cycle_'+cycle+'/'
config.appendDataSearchDir(datadir)

print("\n======= "+inst+" white van reduction =======")

# =================load white van and background and subtract=====================
print("WHITE_VAN "+inst+": Loading white vanadium run# " + str(whitevan))
wv = Load(str(whitevan),LoadMonitors='Exclude')
wv = NormaliseByCurrent('wv')
wv_corrected = Scale('wv',1/whitevan_trans,"Multiply")
if (whitevan_bg is not None):
    print("... subtracting white vanadium background run# %i" % whitevan_bg)
    wv_bg = Load(str(whitevan_bg),LoadMonitors='Exclude')
    wv_bg = NormaliseByCurrent('wv_bg')
    wv_corrected = wv/whitevan_trans - wv_bg
    
# ==============create wv_norm workspace (white vanadium integrals)===============
print("... Calculating white vanadium integrals")
print("... Summing between %.1f < \u03BB < %.1f \u212B" % (wv_lrange[0],wv_lrange[1]))
WV_normalised_integrals = ConvertUnits(wv_corrected,'Wavelength')
WV_normalised_integrals = Rebin(WV_normalised_integrals,str(wv_lrange[0])+',100,'+str(wv_lrange[1]))
wv_normt = Transpose(WV_normalised_integrals); wv_normt = CropWorkspace(wv_normt,XMin=wv_detrange[0],Xmax=wv_detrange[1])
wv_scale = Integration(wv_normt); scale_factor = len(wv_normt.dataY(0))/wv_scale.dataY(0)[0]
WV_normalised_integrals = Scale(WV_normalised_integrals, scale_factor, 'Multiply')

# ===========================output integrals file================================
ofile = 'WV_'+str(whitevan)+'.txt'
print("WHITE_VAN "+inst+": Saving integrals in %s" % ofile)
print("... average intensity = %.2f" % (1/scale_factor))
SaveAscii(WV_normalised_integrals,ofile)

# ===================================cleanup======================================
ws_list = ADS.getObjectNames()
if (not idebug):
    nx_list = [ss for ss in ws_list if 'wv' in ss]
    for ss in nx_list:
        ADS.remove(ss)

# ==================================================================
print("\nWHITE_VAN "+inst+": Reduction complete in %.1f seconds" % (time.time() - t))

