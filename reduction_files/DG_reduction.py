# =======================================================
#
#      EXCITATIONS INSTRUMENTS REDUCTION SCRIPT
#      JRS 3/7/23
#
#      Reads sample run(s), corrects for backgound,
#      normalises to a white vanadium run, and  the
#      absolute normalisation factor for each Ei (if MV
#      file is specified).
#      - Masks bad detectors with hard mask only.
#      - Converts time-of-flight to energy transfer
#      - performs Q-rebinning for QENS data ('_red.nxs' format)
#      - Outputs .nxspe, .nxs or _red.nxs fiies
#
#========================================================

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import matplotlib.pyplot as plt
import numpy as np
import time
from importlib import reload

#=======================User Inputs=======================
powder         = False                    #powder or 1to1 map
sumruns        = False                    #set to True to sum sample runs
sample         = 93368                    #sample runs
sample_bg      = None                     #background runs
wv_file        = 'WV_91329.txt'           #white vanadium integral file (mandatory)
Ei_list        = [22.78,7.52,3.7,2.2]     #incident energies (Ei)
Erange         = [-0.8,0.0025,0.8]        #energy transfer range to output in fractions of Ei
trans          = [1,1,1,1]                #elastic line transmission factors for each Ei
mask           = "LET_mask_231.xml"       #hard mask
#========================================================

#==================Absolute Units Inputs=================
mv_file       = None        # pre-processed MV calibration
monovan_mass  = 3.56        # mass of LET vanadium cylinder
sample_mass   = 1           # mass of the sample
sample_fwt    = 50.9415     # formula weight of sample
#========================================================

#==================Local Contact Inputs==================
# MARI monitors:    m2spec=2,     m3spec=3     - fixei = False
# MERLIN monitors:  m2spec=69636, m3spec=69640 - fixei = False
# MAPS monitors:    m2spec=41475, m3spec=41476 - fixei = False
# LET monitors:     m2spec=98310, m3spec= None - fixei = True

config['default.instrument'] = 'LET'
cycle   = '23_2'                #cycle number
m2spec  = 98310                 #specID of monitor2 (pre-sample)
m3spec  = None                  #specID of monitor3 (post-sample)
fixei   = True                  #True for LET since no monitor 3
powdermap = 'LET_rings_222.xml' #rings mapping file - must be .xml format
file_wait = 30                  #wait for data file to appear (seconds)
keepworkspaces = True           #should be false for Horace scans
saveformat = '.nxspe'           #format of output, ".nxspe", ".nxs"
QENS = False                    #output Q-binned data for QENS data analysis "_red.nxs"
Qbins = 20                      #approximate number of Q-bins (QENS)
theta_range = [5.,65.]          #useful theta range for Q-binning (QENS)

idebug  = False                 #keep workspaces and check absolute units on elastic line
#========================================================

#========================================================
def tryload(irun):             #loops till data file appears
    while True:
        try:
            ws = Load(str(irun),LoadMonitors=True)
        except:
            print("...waiting for run #%i" % irun)
            time.sleep(file_wait)
            continue
        break
#========================================================

# ==============================setup directroties================================
inst = config['default.instrument']
if inst == 'MARI':
    source = 'Moderator'
else:
    source = 'undulator'

datadir = '/archive/NDX'+inst+'/instrument/data/cycle_'+cycle+'/'
mapdir  = '/usr/local/mprogs/InstrumentFiles/'+inst.swapcase()+'/'
config.appendDataSearchDir(mapdir)
config.appendDataSearchDir(datadir)
ws_list = ADS.getObjectNames()

print("\n======= "+inst+" data reduction =======")
print('Working directory... %s\n' % ConfigService.Instance().getString('defaultsave.directory'))

# ============================create lists if necessary==========================
if type(sample) is int:
    sample = [sample]
if sample_bg is not None and type(sample_bg) is int:
    sample_bg = [sample_bg]

#==================================load hard mask================================
if mask is None:
    print(inst+": WARNING - No hard mask!  Bad detectors wont be removed")
if mask not in ws_list and mask is not None:
    print(inst+": Loading hard mask - %s" % mask)
    LoadMask(inst,mask,OutputWorkspace=mask)
else:
    print(inst+": Using previously loaded hard mask - %s" % mask)

# ===============================load whitevan file==============================
if wv_file is None:
    print(inst+": ERROR - white vanadium calibration file missing")
    exit(1)
if wv_file not in ws_list:
    print(inst+": Loading white vanadium - %s" % wv_file)
    LoadAscii(Filename=wv_file,OutputWorkspace=wv_file)
    ReplaceSpecialValues(wv_file, SmallNumberThreshold=1e-20, SmallNumberValue='NaN',OutputWorkspace=wv_file)
else:
    print(inst+": Using previously loaded white vanadium - %s" % wv_file)

# ===============================load monovan file===============================
mv_fac = []
if mv_file is not None:
    if mv_file not in ws_list:
        print(inst+": Loading monovan calibration factors - %s" % mv_file)
        LoadAscii(mv_file,OutputWorkspace=mv_file)
    mv_ws = mtd[mv_file]
    mv_eis = mv_ws.readX(0)
    mv_cal = mv_ws.readY(0)
# check that monovan is compatible with Ei_list)
    Ei_diff = sum([x-y for x,y in zip(mv_eis,Ei_list)])
    if (abs(Ei_diff) > 0):
        print("----ERROR: Monovan file Eis not compatible with Ei_list")
    for Ei in Ei_list:
        ii = np.where(mv_eis == Ei)
        mvf = mv_cal[ii][0]
        van_fwt = 50.9415               # vanadium molar mass (g)
        van_xs = 5080. / 4. / np.pi     # vanadium differential cross-section (mbarns/st)
        mvf = (monovan_mass / van_fwt) / (sample_mass / sample_fwt) * van_xs / mvf
        mv_fac.append(mvf)
else:
    print(inst+": Skipping absolute calibration")
    mv_fac = [x/x for x in Ei_list]     # monovan factors = 1 by default

# =======================load background runs and sum=========================
if sample_bg is not None:
    for irun in sample_bg:
        ws_bg = Load(str(irun),LoadMonitors=False)
        if (irun == sample_bg[0]):
            print("background run #%i loaded" % irun)
            w_buf = CloneWorkspace('ws_bg')
        else:
            print("background run #%i added" % irun)
            w_buf = Plus('w_buf', 'ws_bg')
    ws_bg = CloneWorkspace('w_buf')
    ws_bg = NormaliseByCurrent('ws_bg')

# =======================sum sample runs if required=========================
if sumruns:
    for irun in sample:
        tryload(irun)
        if (irun == sample[0]):
            print("run #%i loaded" % irun)
            w_buf = CloneWorkspace('ws'); w_buf_monitors = CloneWorkspace('ws_monitors')
        else:
            print("run #%i added" % irun)
            w_buf = Plus('w_buf', 'ws'); w_buf_monitors = Plus('w_buf_monitors', 'ws_monitors')
    ws = CloneWorkspace('w_buf'); ws_monitors = CloneWorkspace('w_buf_monitors')
    sample = [sample[0]]

# =======================sample run loop=====================================
for irun in sample:
    t = time.time()         #start the clock

#check loaded workspaces and remove old nxspe workspaces
    if not keepworkspaces:
        ws_list = ADS.getObjectNames()
        nx_list = [ss for ss in ws_list if saveformat in ss]
        for ss in nx_list:
            ADS.remove(ss)

    print("============")
    if not sumruns:
        tryload(irun)
        print("Loading run# %i" % irun)
    if inst == 'MARI':
        ws = RemoveSpectra('ws',[0])
    ws = NormaliseByCurrent('ws')

# ============================= Ei loop =====================================
    for ienergy in range(len(Ei_list)):
       Ei  = Ei_list[ienergy]
        origEi = Ei
        tr  = trans[ienergy]
        mvf = mv_fac[ienergy]
        print("\n"+inst+": Reducing data for Ei=%.2f meV" % Ei)

        ws_corrected = Scale('ws',1/tr,"Multiply")
        if sample_bg is not None:
            print("... subtracting background - transmission factor = %.2f " % tr)
            ws_corrected  = ws/tr - ws_bg

# normalise to WB vanadium and apply fixed mask
        print("... normalising/masking data")
        ws_norm = Divide('ws_corrected',wv_file)          # white beam normalisation
        MaskDetectors(ws_norm,MaskedWorkspace=mask,ForceInstrumentMasking=True)

# t2e section
        print("... t2e section")
        ws_monitors = mtd['ws_monitors']
        spectra = ws_monitors.getSpectrumNumbers()
        index = spectra.index(m2spec)
        m2pos = ws.detectorInfo().position(index)[2]

# this section shifts the time-of-flight such that the monitor2 peak
# in the current monitor workspace (post monochromator) is at t=0 and L=0
# note that the offest is not the predicted value due to energy dependence of the source position

        if m3spec is not None and not fixei:
            (Ei,mon2_peak,_,_) = GetEi(ws_monitors,Monitor1Spec=m2spec,Monitor2Spec=m3spec,EnergyEstimate=Ei)
            print("... refined Ei=%.2f meV" % Ei)
        else:
            (Ei,mon2_peak,_,_) = GetEi(ws_monitors,Monitor2Spec=m2spec,EnergyEstimate=Ei,FixEi=fixei)

        print("... m2 tof=%.2f mus, m2 pos=%.2f m" % (mon2_peak,m2pos))

        ws_norm = ScaleX(ws_norm, Factor=-mon2_peak, Operation='Add', InstrumentParameter='DelayTime', Combine=True)
        MoveInstrumentComponent(ws_norm, ComponentName=source, Z=m2pos, RelativePosition=False)

        ws_out = ConvertUnits(ws_norm,'DeltaE',EMode='Direct',EFixed=Ei)
        ws_out = Rebin(ws_out, [x*Ei for x in Erange], PreserveEvents=False)
        ws_out = DetectorEfficiencyCor(ws_out, IncidentEnergy=Ei)
        ws_out = CorrectKiKf(ws_out, Efixed=Ei, EMode='Direct')

# monovan scaling
        if mv_file is not None:
            print("... applying mono van calibration factor %.1f " % mvf)
        ws_out = Scale('ws_out',mvf,'Multiply')

# rings grouping if desired
        ofile_suffix='_1to1'
        if powder or inst == 'MARI' or QENS:
            ws_out=GroupDetectors(ws_out, MapFile=powdermap, Behaviour='Average')
            ofile_suffix = '_powder'
            if inst == 'MARI' or QENS:
                ofile_suffix = ''
            print("... powder grouping using %s" % powdermap)

# output nxspe file
        ofile = '{:s}{:d}_{:g}meV{:s}'.format(inst[:3],irun,origEi,ofile_suffix)

# check elastic line (debug mode)
        if idebug:
            Rebin('ws_out',[-0.05*Ei,100,0.05*Ei], PreserveEvents=False, OutputWorkspace=ofile+'_elastic')
            Transpose(ofile+'_elastic',OutputWorkspace=ofile+'_elastic')

# output appropriate formats
        ws_out = ConvertToDistribution(ws_out)
        if QENS:
            saveformat = '.nxs'
        print(inst+": Writing %s" % ofile+saveformat)
        if saveformat.lower() == '.nxspe':
            SaveNXSPE('ws_out',ofile+saveformat,Efixed=Ei,KiOverKfScaling=True)
        elif saveformat.lower() == '.nxs':
            SaveNexus('ws_out', ofile+saveformat)
        if QENS:
            print("... outputting QENS '_red' format")
            theta = np.array([theta_range[0],theta_range[1]])*np.pi/180.
            Q     = 1.39 * np.sqrt(Ei) * np.sin(theta)
            Q     = np.around(Q*10) / 10.
            Qbin  = int((Q[1] - Q[0])) / Qbins
            print("... Q-axis = [%g,%g,%g]" % (Q[0]+Qbin,Qbin,Q[1]-Qbin))
            ws_out = SofQW3('ws_out', [Q[0]+Qbin,Qbin,Q[1]-Qbin], "Direct", Efixed=Ei)
# these lines are Anthony Lim's method of removing NaNs
# NaN are changed to zeros, and then placed at the end of the energy range
            spectra = range(ws_out.getNumberHistograms())
            for ispec in spectra:
                x = ws_out.dataX(ispec)
                y = ws_out.dataY(ispec)
                e = ws_out.dataE(ispec)
                smidge = 1e-6
                for iq in range(len(x)-1):
                    if np.isnan(y[iq]):
                        x[iq]=np.max(x) + smidge
                        y[iq]=0.0
                        e[iq]=0.0
                        smidge += 1e-6
            SaveNexus('ws_out', ofile+"_red"+saveformat)

        CloneWorkspace('ws_out',OutputWorkspace=ofile+saveformat)

# ============================= End of Ei loop ================================

    print("\n"+inst+": Reduction complete in %.1f seconds\n" % (time.time() - t))

# ============================= End of run loop ================================

# cleanup
if not idebug:
    ws_list = ADS.getObjectNames()
    nx_list = [ss for ss in ws_list if 'w_buf' in ss or 'ws' in ss]
    for ss in nx_list:
        ADS.remove(ss)
