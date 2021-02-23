from Bio import AlignIO
import numpy as np
import os


class MSA:

    """
    Receive a MSA .fas alignment, parse and return the positions of information channels array.

    """

    def __init__(self):

        self.__MSA_A_file = None
        self.__MSA_B_file = None
        self.__MSA_A = None
        self.__MSA_B = None
        self.__encoded_MSA_A = np.zeros((0, 0))
        self.__encoded_MSA_B = np.zeros((0, 0))
        self.__aa_code = {"A": 0, "R": 1, "N": 2, "D": 3, "Q": 4,
                          "E": 5, "G": 6, "H": 7, "L": 8, "K": 9,
                          "M": 10, "F": 11, "S": 12, "T": 13, "W": 14,
                          "Y": 15, "C": 16, "I": 17, "P": 18, "V": 19,
                          "-": 20, ".": 20, "B": 2, "Z": 4, "X": 20, "J": 20}
        self.__q = 21
        self.__sequence_weight = np.zeros([0, 0])

    def __set_encoded(self, msa):

        encoded = np.empty((len(msa), len(msa[0])), dtype=np.int16)

        for idx, aa in np.ndenumerate(msa):
            aa = aa.upper()
            encoded[idx] = int(self.__aa_code[aa])

        return encoded

    def set_MSA_A(self, file_name):

        try:
            self.__MSA_A = AlignIO.read(file_name, "fasta")
        except:
            return False

        self.__MSA_A_file = str(file_name)
        if self.__MSA_B and len(self.__MSA_B) != len(self.__MSA_A):
            return False
        else:
            return True

    def set_MSA_B(self, file_name):
        
        try:
            self.__MSA_B = AlignIO.read(file_name, "fasta")
        except:
            return False

        self.__MSA_B_file = str(file_name)
        if self.__MSA_A and len(self.__MSA_A) != len(self.__MSA_B):
            return False
        else:
            return True

    def sub_full_MSA(self, ic_list, swap=None):

        if not self.__encoded_MSA_A.size:
            self.__encoded_MSA_A = self.__set_encoded(self.__MSA_A)
        if not self.__encoded_MSA_B.size:
            self.__encoded_MSA_B = self.__set_encoded(self.__MSA_B)

        if swap is not None:
            encoded_b = self.__encoded_MSA_B[swap.flatten()]
        else:
            encoded_b = self.__encoded_MSA_B

        encoded_msa = np.take(np.concatenate((self.__encoded_MSA_A, encoded_b), axis=1), ic_list, axis=1)
        return encoded_msa

    @property
    def ref_msa_A(self):

        a = self.__MSA_A[0].format("fasta").rsplit()[1:]
        a = "".join(a)
        return a

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

    @property
    def MSA_A(self):

        if not self.__encoded_MSA_A.size:
            self.__encoded_MSA_A = self.__set_encoded(self.__MSA_A)
        return self.__encoded_MSA_A

    def sub_MSA_A(self, ICs):

        if not self.__encoded_MSA_A.size:
            self.__encoded_MSA_A = self.__set_encoded(self.__MSA_A)
        sub = np.take(self.__encoded_MSA_A, ICs, axis=1)
        return sub

    @property
    def MSA_A_sz(self):
        return self.__encoded_MSA_A.shape[1]

    @property
    def ref_msa_B(self):

        b = self.__MSA_B[0].format("fasta").rsplit()[1:]
        b = "".join(b)
        return b

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
    def MSA_B(self):

        if not self.__encoded_MSA_B.size:
            self.__encoded_MSA_B = self.__set_encoded(self.__MSA_B)
        return self.__encoded_MSA_B

    def sub_MSA_B(self, ICs, swap=None):

        if not self.__encoded_MSA_B.size:
            self.__encoded_MSA_B = self.__set_encoded(self.__MSA_B)
        sub = np.take(self.__encoded_MSA_B, ICs, axis=1)

        if swap:
            sub[swap.flatten()] = sub[np.arange(self.__encoded_MSA_B.shape[0])]

        return sub

    @property
    def MSA_B_sz(self):
        return self.__encoded_MSA_B.shape[1]

    @property
    def num_of_seqs(self):
        return self.__encoded_MSA_A.shape[0]

    @property
    def full_encoded(self):

        if not self.__encoded_MSA_A.size:
            self.__encoded_MSA_A = self.__set_encoded(self.__MSA_A)
        if not self.__encoded_MSA_B.size:
            self.__encoded_MSA_B = self.__set_encoded(self.__MSA_B)
        return np.concatenate((self.__encoded_MSA_A, self.__encoded_MSA_B), axis=1)

    @property
    def sequence_weight(self):
        return self.__sequence_weight

    @sequence_weight.setter
    def sequence_weight(self, vector):
        self.__sequence_weight = np.copy(vector)

    @property
    def Meff(self):
        return np.sum(self.sequence_weight)

    @property
    def q(self):
        return self.__q
