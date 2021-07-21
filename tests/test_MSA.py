import unittest
import numpy as np
import sys 
sys.path.append('../coevtools')
from Bio import AlignIO
from coevtools.MSA import MSA
from AB_1BXR.PDB_std import A_ICs_AB, B_ICs_AB, B_ICs_shift


class MSATestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("Initializing MSA tests:")

    def setUp(self):
        self.msa = MSA("AB_1BXR/MSA_std/1BXR_A.fas", "AB_1BXR/MSA_std/1BXR_B.fas")
        print("\n{} ".format(self.shortDescription()))

    def test_encoded_MSA_A(self):
        """ encoding MSA A """
        std = np.load("AB_1BXR/MSA_std/A_encoded.npy")
        self.assertTrue(np.array_equal(self.msa.MSA_A(), std))

    def test_encoded_MSA_B(self):
        """ encoding MSA B """
        self.subTest("encoded without swapping")
        std = np.load("AB_1BXR/MSA_std/B_encoded.npy")
        self.assertTrue(np.array_equal(self.msa.MSA_B(), std))

    def test_full_encoded_AB(self):
        """ encoding full MSA """
        std = np.load("AB_1BXR/MSA_std/full_encoded.npy")
        self.assertTrue(np.array_equal(self.msa.MSA_AB(), std))

    def test_sub_MSA_A(self):
        """ extracting sub MSA A """
        A_ICs = A_ICs_AB.A_ICs
        encoded = self.msa.MSA_A(ICs=A_ICs)
        std = np.load("AB_1BXR/MSA_std/A_sub_encoded.npy")
        self.assertTrue(np.array_equal(encoded, std))

    def test_sub_MSA_B(self):
        """ extracting sub MSA B """
        B_ICs = B_ICs_AB.B_ICs
        encoded = self.msa.MSA_B(ICs=B_ICs)
        std = np.load("AB_1BXR/MSA_std/B_sub_encoded.npy")
        self.assertTrue(np.array_equal(encoded, std))

    def test_sub_full_encoded_AB(self):
        """ extracting full sub MSA """
        A_ICs = A_ICs_AB.A_ICs
        B_ICs = B_ICs_shift.B_ICs_shift
        with self.subTest("encoded without swapping"):
            encoded = self.msa.MSA_AB(ICs=A_ICs + B_ICs)
            std = np.load("AB_1BXR/MSA_std/sub_full_encoded.npy")
            self.assertTrue(np.array_equal(encoded, std))
        with self.subTest("encoded with swapping"):
            idx = np.load("AB_1BXR/MSA_std/vec_swap_row_1.npy")
            std = np.load("AB_1BXR/MSA_std/arr_swap_B_row_vec.npy")
            encoded = self.msa.MSA_AB(ICs=A_ICs + B_ICs, swap=idx)
            self.assertTrue(np.array_equal(encoded, std))

    def test_msa_A_fasta(self):
        """ extracting fasta from sub MSA A """
        fasta = self.msa.MSA_A_fasta(A_ICs_AB.A_ICs, "test")
        std_fasta = AlignIO.read("AB_1BXR/MSA_std/sub_msa_A_AB.fas", "fasta")
        flag = True
        for i in range(len(std_fasta)):
            if std_fasta[i].seq != fasta[i].seq: flag = False
        self.assertTrue(flag)
    
    def test_msa_B_fasta(self):
        """ extracting fasta from sub MSA B """
        fasta = self.msa.MSA_B_fasta(B_ICs_AB.B_ICs, "test")
        std_fasta = AlignIO.read("AB_1BXR/MSA_std/sub_msa_B_AB.fas", "fasta")
        flag = True
        for i in range(len(std_fasta)):
            if std_fasta[i].seq != fasta[i].seq: flag = False
        self.assertTrue(flag)

    def test_sequence_weight(self):
        """ calculating sequence weight - theta = 0.8 """
        seq_weight = self.msa.calc_sequence_weight(theta=0.8)
        std = np.load("AB_1BXR/MSA_std/l05t08_seqweight.npy")
        self.assertTrue(np.array_equal(seq_weight, std))


def suite():
    suite = unittest.TestSuite()
    suite.addTest(MSATestCase('test msa input'))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
