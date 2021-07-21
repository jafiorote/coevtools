import unittest
import sys
sys.path.append('../coevtools')
import numpy as np
from coevtools.PDB import PDB
from coevtools.MSA import MSA
from coevtools.Information import Information

class CoevolutionTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("Initializing Coevolution functions tests:")

    def setUp(self):
        self.pdb = PDB("AB_1BXR/pdbs/1BXR_AB.pdb")
        self.msa = MSA("AB_1BXR/MSA_std/1BXR_A.fas", "AB_1BXR/MSA_std/1BXR_B.fas", theta=0.8)
        print("\n{} ".format(self.shortDescription()))

    def test_site_freq(self):
        """ Calculating site frequencies """
        # standard lambda  = 0.5
        ics = self.pdb.ICs
        sf = Information.site_freq(self.msa, ics)
        std = np.load("AB_1BXR/Information_std/l05t08_sitefreq.npy")
        self.assertTrue(np.array_equal(std, sf))

    def test_pair_freq(self):
        """ Calculating pair frequencies """
        # standard lambda  = 0.5
        ics = self.pdb.ICs
        std = np.load("AB_1BXR/Information_std/l05t08_pairfreq.npy")
        sf = Information.site_freq(self.msa, ics)
        pf = Information.pair_freq(self.msa, ics, sf)
        self.assertTrue(np.array_equal(std, pf))

    def test_shannon_information(self):
        """ Calculating Shannon Information"""
        std_mi = np.load("AB_1BXR/Information_std/l05t08_mi_matrix.npy")
        std_h = np.load("AB_1BXR/Information_std/l05t08_h_matrix.npy")

        sf = Information.site_freq(self.msa, self.pdb.ICs)
        pf = Information.pair_freq(self.msa, self.pdb.ICs, sf)
        sh = Information.shannon_info(self.msa, self.pdb, sf, pf)
        self.assertTrue(np.array_equal(std_mi, sh[0]))
        self.assertTrue(np.array_equal(std_h, sh[1]))
        #self.assertAlmostEqual(np.sum(sh[0]), 0.072239366200623822)

    def test_r_value(self):
        """ Calculating linear correlation coeficient"""
        r = Information.mirror_tree(self.msa.MSA_A_fasta(self.pdb.A_ICs, "teste"),
                                    self.msa.MSA_B_fasta(self.pdb.B_ICs, "teste"))
        self.assertAlmostEqual(r, 0.68620817680883106)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(CoevolutionTestCase('test msa input'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
