from Bio import AlignIO
import numpy as np
from scipy import spatial
import os
import sys


class MSA:

    """
    Receive a MSA .fas alignment, parse and return the positions of information channels.

    """

    def __init__(self, MSA_A_file, MSA_B_file, theta=1.0):

        self.__MSA_A_file = MSA_A_file
        self.__MSA_B_file = MSA_B_file
        self.__MSA_A = []
        self.__MSA_B = []
        self.__msa_consistency()
        self.__aa_code = {"A": 0, "R": 1, "N": 2, "D": 3, "Q": 4,
                          "E": 5, "G": 6, "H": 7, "L": 8, "K": 9,
                          "M": 10, "F": 11, "S": 12, "T": 13, "W": 14,
                          "Y": 15, "C": 16, "I": 17, "P": 18, "V": 19,
                          "-": 20, ".": 20, "B": 2, "Z": 4, "X": 20, "J": 20}
        self.__encoded_MSA_A = self.__set_encoded(self.__MSA_A)
        self.__encoded_MSA_B = self.__set_encoded(self.__MSA_B)
        self.__q = 21
        self.__sequence_weight = self.calc_sequence_weight(theta)

    def __set_encoded(self, msa):

        """
        Encoded MSA into a numpy array.

        """

        encoded = np.empty((len(msa), len(msa[0])), dtype=np.int16)

        for idx, aa in np.ndenumerate(msa):
            aa = aa.upper()
            encoded[idx] = int(self.__aa_code[aa])

        return encoded

    def __msa_consistency(self):

        """
        Check MSA files.

        """

        try:
            self.__MSA_A = AlignIO.read(self.__MSA_A_file, "fasta")
        except OSError:
            print(f"Could not open/read file: {self.__MSA_A_file}")
            sys.exit()

        try:
            self.__MSA_B = AlignIO.read(self.__MSA_B_file, "fasta")
        except OSError:
            print(f"Could not open/read file: {self.__MSA_B_file}")
            sys.exit()

        if len(self.__MSA_A) != len(self.__MSA_B):
            print("MSA A and MSA B must have same sequences number.")
            sys.exit()

    def calc_sequence_weight(self, theta):

        """
        Calculate the weight of individual seqs in MSA.

        """

        hamming_distance = spatial.distance.pdist(self.MSA_AB(), "hamming")
        weight_matrix = spatial.distance.squareform(hamming_distance < (1.0 - theta))
        return 1.0 / (np.sum(weight_matrix, axis=1) + 1.0)

    def MSA_AB(self, ICs=None, swap=None):

        """
        Return encoded MSA.
        ____________
        parameters:
        ICs: Information channels. Columns to take from full MSA.
        swap:

        """

        encoded_b = self.__encoded_MSA_B

        if swap is not None:
            encoded_b = self.__encoded_MSA_B[swap.flatten()]
        if ICs is not None:
            encoded_msa = np.take(np.concatenate((self.__encoded_MSA_A, encoded_b), axis=1), ICs, axis=1)
        else:
            encoded_msa = np.concatenate((self.__encoded_MSA_A, encoded_b), axis=1)

        return encoded_msa

    def MSA_A(self, ICs=None, swap=None):

        msa_a = self.__encoded_MSA_A

        if ICs:
            msa_a = np.take(self.__encoded_MSA_A, ICs, axis=1)
        if swap:
            msa_a[swap.flatten()] = msa_a[np.arange(self.__encoded_MSA_A.shape[0])]

        return msa_a

    def MSA_B(self, ICs=None, swap=None):

        msa_b = self.__encoded_MSA_B

        if ICs:
            msa_b = np.take(self.__encoded_MSA_B, ICs, axis=1)
        if swap:
            msa_b[swap.flatten()] = msa_b[np.arange(self.__encoded_MSA_B.shape[0])]

        return msa_b

    def MSA_A_fasta(self, ics, name):

        temp = ""

        for seq in self.__MSA_A:
            temp += ">" + seq.id + "\n"
            line = ""
            for i in ics:
                line += seq[i]
            temp += line + "\n"
        with open("sub_msa_A_{}.fas".format(name), 'w') as fl:
            fl.write(temp)
        fl.close()

        msa_A = AlignIO.read("sub_msa_A_{}.fas".format(name), "fasta")
        os.remove("sub_msa_A_{}.fas".format(name))
        return msa_A

    def MSA_B_fasta(self, ics, name):

        temp = ""

        for seq in self.__MSA_B:
            temp += ">" + seq.id + "\n"
            line = ""
            for i in ics:
                line += seq[i]
            temp += line + "\n"
        with open("sub_msa_B_{}.fas".format(name), 'w') as fl:
            fl.write(temp)
        fl.close()

        msa_B = AlignIO.read("sub_msa_B_{}.fas".format(name), "fasta")
        os.remove("sub_msa_B_{}.fas".format(name))
        return msa_B

    @property
    def MSA_A_sz(self):
        return self.__encoded_MSA_A.shape[1]

    @property
    def MSA_B_sz(self):
        return self.__encoded_MSA_B.shape[1]

    @property
    def num_of_seqs(self):
        return self.__encoded_MSA_A.shape[0]

    @property
    def sequence_weight(self):
        return self.__sequence_weight

    @property
    def Meff(self):
        return np.sum(self.sequence_weight)

    @property
    def q(self):
        return self.__q
