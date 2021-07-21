import unittest
import numpy as np
import sys
sys.path.append('../coevtools')
from coevtools.PDB import PDB
from AB_1BXR.PDB_std import A_ICs_AB, B_ICs_AB, B_ICs_shift


class PDBTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("Initializing Model tests:")

    def setUp(self):
        self.pdb = PDB("AB_1BXR/pdbs/1BXR_AB.pdb")
        print("\n{} ".format(self.shortDescription()))

    def test_contact_map(self):
        """ localizating interface residues """
        std = np.sort(np.load("AB_1BXR/PDB_std/c_map.npy"), axis=None)
        c_map = np.sort(self.pdb.pairs_shift, axis=None)
        self.assertTrue(np.array_equal(c_map, std))
        
    def test_ICs(self):
        """ selecting information channels """
        with self.subTest():
            self.assertEqual(set(self.pdb.A_ICs), set(A_ICs_AB.A_ICs))
        with self.subTest():
            self.assertEqual(set(self.pdb.B_ICs), set(B_ICs_AB.B_ICs))
        with self.subTest():
            self.assertEqual(set(self.pdb.shift_B_ICs), set(B_ICs_shift.B_ICs_shift))
        with self.subTest():
            self.assertEqual(set(self.pdb.ICs), set(A_ICs_AB.A_ICs + B_ICs_shift.B_ICs_shift))

def suite():
    suite = unittest.TestSuite()
    suite.addTest(PDBTestCase('test model parsing'))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
