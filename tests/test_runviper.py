## Licensed under a GPLv3 style license - see LICENSE

from astropy.io import fits
import pytest
import unittest
import importlib
import numpy as np
import os
import subprocess
from cpl.core import Table

path = os.getcwd()
directory = os.path.dirname(os.path.realpath(__file__)) + os.sep 

Inst = importlib.import_module('inst.inst_CRIRES')
Tpl = Inst.Tpl

try:
    # remove existing files from last test run
    os.remove(directory + "tmp_rvo_par.fits")
    os.remove(directory + "tmp1_rvo_par.fits")
    os.remove(directory + "tmp_tpl.fits")
except:
    pass


#@pytest.fixture(scope="class")
#def run_viper(args1):
#    return os.system(args1)


class test_viper(unittest.TestCase):
    def test_pycpl(self):
        # Test if PyCPL is installed (just CRIRES+)
        # viper will make use of astropy if PyCPL is not available
        
        self.assertEqual(Inst.pycpl, 1)

    def test_template(self):
        # Test run to create a telluric free template out of several observations
        # for testing just two orders and two observations were used

        os.system("python3 -m viper " + directory + "'test_data/22*' -inst CRIRES -deg_norm 0 -deg_wave 2 -deg_bkg 1 -oversampling 1 -createtpl -telluric add -tsig 10 -nocell -nset :2 -oset 7,12 -output_format cpl -tag " + directory + "tmp")

        # test if template was generated
        assert os.path.exists(directory + "tmp_tpl.fits")

        # test if results fit the test_template
        wave_tpl0, spec_tpl0 = Tpl(directory + "test_compare/test_tpl.fits", order=7)
        wave_tpl, spec_tpl = Tpl(directory + "tmp_tpl.fits", order=7)
        assert np.isclose(wave_tpl0[1000], wave_tpl[1000], rtol=1e-3, atol=0)
        assert np.isclose(spec_tpl0[1000], spec_tpl[1000], rtol=1e-3, atol=0)

        wave_tpl0, spec_tpl0 = Tpl(directory + "test_compare/test_tpl.fits", order=12)
        wave_tpl, spec_tpl = Tpl(directory + "tmp_tpl.fits", order=12)
        assert np.isclose(wave_tpl0[1000], wave_tpl[1000], rtol=1e-3, atol=0)
        assert np.isclose(spec_tpl0[1000], spec_tpl[1000], rtol=1e-3, atol=0)

    def test_calc_RV(self):
        # Calculte RVs with test template
        # just two orders of two observations for fast test

        os.system("python3 -m viper " + directory + "'test_data/SGC*' " + directory + "test_compare/test_tpl.fits -inst CRIRES -deg_norm 2 -deg_wave 2 -deg_bkg 1 -tsig 1 -telluric add -tellshift -kapsig 0 4.5 -oversampling 1 -nset :2 -oset 7,12 -output_format cpl -tag " + directory + "tmp1")

        # test if all products were generated
        assert os.path.exists(directory+"tmp1_rvo_par.fits")

        # test reults; comparison to pre-generated RV results
        tbl = Table.load(directory+"test_compare/test_rvo_par.fits", 1)
        print('----',  np.array(tbl["BJD"])[0])
        bjd, RV, e_RV = np.array(tbl["BJD"])[0], np.array(tbl["RV"])[0], np.array(tbl["e_RV"])[0]

        tbl = Table.load(directory+"tmp1_rvo_par.fits", 1)
        bjd2, RV2, e_RV2 = np.array(tbl["BJD"])[0], np.array(tbl["RV"])[0], np.array(tbl["e_RV"])[0]

        assert np.isclose(bjd, bjd2, rtol=1e-5, atol=0)
        assert np.isclose(RV, RV2, rtol=1e-1, atol=0)
        assert np.isclose(e_RV, e_RV2, rtol=1e-1, atol=0)


def main():
    unittest.main()

if __name__ == "__main__":
    main()
