#!/usr/bin/env python3
import unittest
import numpy as np
import os, sys, re
import mantid.simpleapi as s_api
from os.path import abspath, dirname
from mantid.simpleapi import LoadNexusProcessed

class DGReductionTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Sets up Mantid paths
        cls.repopath = abspath(os.path.join(dirname(__file__), '..'))
        cls.scriptpath = os.path.join(cls.repopath, 'reduction_files')
        cls.datapath = os.path.join(cls.repopath, 'datafiles')
        cls.instpath = os.path.join(cls.repopath, 'InstrumentFiles')
        cls.outputpath = os.path.join(cls.repopath, 'tests')
        sys.path.append(cls.outputpath)
        sys.path.append(cls.scriptpath)
        s_api.config['defaultsave.directory'] = cls.outputpath
        s_api.config.appendDataSearchDir(cls.datapath)
        for inst in ['let', 'maps', 'mari', 'merlin']:
            s_api.config.appendDataSearchDir(os.path.join(cls.instpath, inst))


    @staticmethod
    def substitute_file(filename, out_file, subs_dict):
        with open(filename, 'r') as f:
            ftxt = f.read()
        for ky, val in subs_dict.items():
            ftxt = re.sub(ky, val, ftxt)
        with open(out_file, 'w') as f:
            f.write(ftxt)


    #@classmethod
    #def tearDownClass(cls):
    #    # Removes all generated files


    def test_MARI(self):
        # Checks that the reduction script runs for MARI and generates output files
        subsdict = {'\nconfig':'\n#config',
                    'save_dir = ':'save_dir = None #',
                    'INSTRUMENT_NAME':'MARI',
                    'MASK_FILE_XML':'mari_mask2023_1.xml',
                    'RINGS_MAP_XML':'mari_res2013.map',
                    'whitevan\s*=\s*[0-9]*':'whitevan = 28580',
                    'sample\s*=\s*\\[*[\\]0-9,]+':'sample = [28581]',
                    'sample_bg\s*=\s*\\[*[\\]0-9,]+':'sample_bg = None',
                    'wv_file\s*=\s*[\\\'A-z0-9\\.]*':'wv_file = \'WV_28580.txt\'',
                    'wv_detrange\s*=\s*[\\[\\]0-9,]*':'wv_detrange = None',
                    'Ei_list\s*=\s*[\\[\\]\\.0-9,]+.*':'Ei_list = [180, 29.8, 11.7]'}
        s_api.config['default.instrument'] = 'MARI'
        infile = os.path.join(self.scriptpath, 'DG_whitevan.py')
        outfile = os.path.join(self.outputpath, 'mari_whitevan.py')
        self.substitute_file(infile, outfile, subsdict)
        import mari_whitevan
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'WV_28580.txt')))
        infile = os.path.join(self.scriptpath, 'DG_reduction.py')
        outfile = os.path.join(self.outputpath, 'mari_reduction.py')
        self.substitute_file(infile, outfile, subsdict)
        import mari_reduction
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'MAR28581_180meV.nxspe')))
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'MAR28581_29.8meV.nxspe')))
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'MAR28581_11.7meV.nxspe')))
        #Load automatically calls `LoadNeXus` first which chokes on the file...
        #ws = s_api.Load(os.path.join(self.outputpath, 'MAR28581_180meV.nxspe'))


    def test_MARI_multiple_lowE(self):
        subsdict = {'\nconfig':'\n#config',
                    'save_dir = ':'save_dir = None #',
                    'INSTRUMENT_NAME':'MARI',
                    'MASK_FILE_XML':'mari_mask2023_1.xml',
                    'RINGS_MAP_XML':'mari_res2013.map',
                    'whitevan\s*=\s*[0-9]*':'whitevan = 28580',
                    'sample\s*=\s*\\[*[\\]0-9,]+':'sample = [28727, 28728]',
                    'sample_bg\s*=\s*\\[*[\\]0-9,]+':'sample_bg = None',
                    'wv_file\s*=\s*[\\\'A-z0-9\\.]*':'wv_file = \'WV_28580.txt\'',
                    'wv_detrange\s*=\s*[\\[\\]0-9,]*':'wv_detrange = None',
                    'Ei_list\s*=\s*[\\[\\]\\.0-9,]+.*':'Ei_list = [1.84, 1.1]'}
        infile = os.path.join(self.scriptpath, 'DG_reduction.py')
        outfile = os.path.join(self.outputpath, 'mari_reduction_lowE.py')
        self.substitute_file(infile, outfile, subsdict)
        import mari_reduction_lowE
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'MAR28727_1.84meV.nxspe')))
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'MAR28727_1.1meV.nxspe')))


    def test_existing_workspace(self):
        subsdict = {'\nconfig':'\n#config',
                    'save_dir = ':'save_dir = None #',
                    'INSTRUMENT_NAME':'MARI',
                    'MASK_FILE_XML':'mari_mask2023_1.xml',
                    'RINGS_MAP_XML':'mari_res2013.map',
                    'whitevan\s*=\s*[0-9]*':'whitevan = 28580',
                    'sample\s*=\s*\\[*[\\]0-9,]+':'sample = ["ws_existing"]',
                    'sample_bg\s*=\s*\\[*[\\]0-9,]+':'sample_bg = None',
                    'wv_file\s*=\s*[\\\'A-z0-9\\.]*':'wv_file = \'WV_28580.txt\'',
                    'wv_detrange\s*=\s*[\\[\\]0-9,]*':'wv_detrange = None',
                    'Ei_list\s*=\s*[\\[\\]\\.0-9,]+.*':'Ei_list = [1.84, 1.1]'}
        ws_existing = s_api.Load('MAR28728.raw', LoadMonitors='Exclude')
        ws_existing = s_api.RemoveSpectra(ws_existing, [0])
        infile = os.path.join(self.scriptpath, 'DG_reduction.py')
        outfile = os.path.join(self.outputpath, 'mari_reduction_existing.py')
        self.substitute_file(infile, outfile, subsdict)
        import mari_reduction_existing
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'MARws_existing_1.84meV.nxspe')))
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'MARws_existing_1.1meV.nxspe')))


    def test_LET_QENS(self):
        # Checks that the reduction script runs for LET QENS and generates output files
        subsdict = {'\nconfig':'\n#config',
                    'save_dir = ':'save_dir = None #',
                    'INSTRUMENT_NAME':'LET',
                    'MASK_FILE_XML':'LET_mask_222.xml',
                    'RINGS_MAP_XML':'LET_rings_222.xml',
                    'whitevan\s*=\s*[0-9]*':'whitevan = 91329',
                    'sample\s*=\s*\\[*[\\]0-9,]+':'sample = 93338',
                    'sample_bg\s*=\s*\\[*[\\]0-9,]+':'sample_bg = 93329',
                    'wv_file\s*=\s*[\\\'A-z0-9\\.]*':'wv_file = \'WV_91329.txt\'',
                    'Ei_list\s*=\s*[\\[\\]\\.0-9,]+.*':'Ei_list = [3.7, 1.77, 1.03]',
                    'QENS = False':'QENS = True'}
        s_api.config['default.instrument'] = 'LET'
        infile = os.path.join(self.scriptpath, 'DG_whitevan.py')
        outfile = os.path.join(self.outputpath, 'let_whitevan.py')
        self.substitute_file(infile, outfile, subsdict)
        import let_whitevan
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'WV_91329.txt')))
        infile = os.path.join(self.scriptpath, 'DG_reduction.py')
        outfile = os.path.join(self.outputpath, 'let_reduction.py')
        self.substitute_file(infile, outfile, subsdict)
        import let_reduction
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'LET93338_3.7meV_red.nxs')))
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'LET93338_1.77meV_red.nxs')))
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'LET93338_1.03meV_red.nxs')))
        ws = s_api.Load(os.path.join(self.outputpath, 'LET93338_3.7meV_red.nxs'))


    def test_LET_same_angle(self):
        subsdict = {'\nconfig':'\n#config',
                    'save_dir = ':'save_dir = None #',
                    'INSTRUMENT_NAME':'LET',
                    'MASK_FILE_XML':'LET_mask_222.xml',
                    'RINGS_MAP_XML':'LET_rings_222.xml',
                    'whitevan\s*=\s*[0-9]*':'whitevan = 91329',
                    'sample\s*=\s*\\[*[\\]0-9,]+':'sample = [92089, 92168]',
                    'sample_bg\s*=\s*\\[*[\\]0-9,]+':'sample_bg = None',
                    'wv_file\s*=\s*[\\\'A-z0-9\\.]*':'wv_file = \'WV_91329.txt\'',
                    'Ei_list\s*=\s*[\\[\\]\\.0-9,]+.*':'Ei_list = [3.7]',
                    'powder\s*=\s*True': 'powder = False',
                    "saveformat = '.nxspe'":"saveformat = '.nxs'",
                    "same_angle_action = 'ignore'":"same_angle_action = 'replace'"}
        s_api.config['default.instrument'] = 'LET'
        infile = os.path.join(self.scriptpath, 'DG_reduction.py')
        outfile = os.path.join(self.outputpath, 'let_reduction_angle.py')
        self.substitute_file(infile, outfile, subsdict)
        import let_reduction_angle
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'LET92089_3.7meV_1to1.nxs')))
        self.assertFalse(os.path.exists(os.path.join(self.outputpath, 'LET92168_3.7meV_1to1.nxs')))
        ws_out = LoadNexusProcessed(os.path.join(self.outputpath, 'LET92089_3.7meV_1to1.nxs'))
        algs = {f'{h.name()}_{ii}':h.getProperties() for ii, h in enumerate(ws_out.getHistory().getAlgorithmHistories())}
        self.assertTrue('Plus' in [v.split('_')[0] for v in algs.keys()])
        loaded_files = [p.value() for al, pp in algs.items() if 'Load' in al for p in pp if p.name() == 'Filename']
        self.assertTrue(any(['92089' in v for v in loaded_files]) and any(['92168' in v for v in loaded_files]))


    def test_MERLIN(self):
        # Checks that the reduction script runs for MERLIN and generates output files
        subsdict = {'\nconfig':'\n#config',
                    'save_dir = ':'save_dir = None #',
                    'INSTRUMENT_NAME':'MERLIN',
                    'MASK_FILE_XML':'mask_211_fix.xml',
                    'RINGS_MAP_XML':'rings_193.map',
                    'whitevan\s*=\s*[0-9]*':'whitevan = 57088',
                    'sample\s*=\s*\\[*[\\]0-9,]+':'sample = 59151',
                    'sample_bg\s*=\s*\\[*[\\]0-9,]+':'sample_bg = None',
                    'wv_file\s*=\s*[\\\'A-z0-9\\.]*':'wv_file = \'WV_57088.txt\'',
                    'Ei_list\s*=\s*[\\[\\]\\.0-9,]+.*':'Ei_list = [150]',
                    'fixei = True':'fixei = False'}
        s_api.config['default.instrument'] = 'MERLIN'
        infile = os.path.join(self.scriptpath, 'DG_whitevan.py')
        outfile = os.path.join(self.outputpath, 'merlin_whitevan.py')
        self.substitute_file(infile, outfile, subsdict)
        import merlin_whitevan
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'WV_57088.txt')))
        infile = os.path.join(self.scriptpath, 'DG_reduction.py')
        outfile = os.path.join(self.outputpath, 'merlin_reduction.py')
        self.substitute_file(infile, outfile, subsdict)
        import merlin_reduction
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'MER59151_150meV_powder.nxspe')))
        ws = s_api.Load(os.path.join(self.outputpath, 'MER59151_150meV_powder.nxspe'))


    def test_MAPS(self):
        # Checks that the reduction script runs for MAPS and generates output files
        subsdict = {'\nconfig':'\n#config',
                    'save_dir = ':'save_dir = None #',
                    'INSTRUMENT_NAME':'MAPS',
                    'MASK_FILE_XML':'4to1_211_msk.xml',
                    'RINGS_MAP_XML':'MAPS_rings.map',
                    'whitevan\s*=\s*[0-9]*':'whitevan = 41272',
                    'sample\s*=\s*\\[*[\\]0-9,]+':'sample = 41335',
                    'sample_bg\s*=\s*\\[*[\\]0-9,]+':'sample_bg = None',
                    'wv_file\s*=\s*[\\\'A-z0-9\\.]*':'wv_file = \'WV_41272.txt\'',
                    'Ei_list\s*=\s*[\\[\\]\\.0-9,]+.*':'Ei_list = [80]',
                    'fixei = True':'fixei = False',
                    'powder\s*=\s*True':'powder = False',
                    'm2spec = 36867':'m2spec = 41475',
                    'm3spec = 36868':'m3spec = 41476'}
        s_api.config['default.instrument'] = 'MAPS'
        infile = os.path.join(self.scriptpath, 'DG_whitevan.py')
        outfile = os.path.join(self.outputpath, 'maps_whitevan.py')
        self.substitute_file(infile, outfile, subsdict)
        import maps_whitevan
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'WV_41272.txt')))
        infile = os.path.join(self.scriptpath, 'DG_reduction.py')
        outfile = os.path.join(self.outputpath, 'maps_reduction.py')
        self.substitute_file(infile, outfile, subsdict)
        import maps_reduction
        self.assertTrue(os.path.exists(os.path.join(self.outputpath, 'MAP41335_80meV_1to1.nxspe')))
        ws = s_api.Load(os.path.join(self.outputpath, 'MAP41335_80meV_1to1.nxspe'))


    def test_auto_iliad(self):
        # Tests that the iliad driver script with Ei='auto' works
        sys.path.append(self.scriptpath)
        wv_file = os.path.join(self.outputpath, 'WV_28580.txt')
        if os.path.exists(wv_file):
            os.remove(wv_file)
        from reduction_utils import iliad
        iliad(28581, ei='auto', wbvan=28580, inst='MARI', hard_mask_file='mari_mask2023_1.xml',
              powdermap='mari_res2013.map')


    def test_func_continuous(self):
        from reduction_utils import run_reduction
        run_reduction(sample=[97138, 97139], Ei_list=[1.7], sumruns=True, wv_file='WV_91329.txt', 
                      inst='LET', mask='LET_mask_222.xml', powdermap='LET_rings_222.xml',
                      cs_block='T_Stick', cs_block_unit='K', cs_bin_size=10)
        for tt in np.arange(197.9, 288, 10):
            self.assertTrue(os.path.exists(os.path.join(self.outputpath, f'LET97138_1.7meV_{tt}K_powder.nxspe')))
        

if __name__ == '__main__':
    unittest.main()
