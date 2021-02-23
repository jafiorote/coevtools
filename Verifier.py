from Bio import AlignIO
import Levenshtein


class Verifier:

    """ 
    Verify and correct the iteration of the PDB chain. 
    
    """

    @classmethod
    def validate_models(cls, pdbs):
       
        for model in pdbs:

            AIsEqual = cls.identity(pdbs[0].fasta[0], model.fasta[0])
            BIsEqual = cls.identity(pdbs[0].fasta[1], model.fasta[1])
            if not AIsEqual or not BIsEqual:
                return False

        return True

    @staticmethod 
    def identity(seq1, seq2, cutoff=0.98):

        """ Calculate distance between sequences. """

        if len(seq1) != len(seq2):
            raise ValueError("Sequences don't have the same size.")

        sim = Levenshtein.ratio(seq1, seq2)

        return True if sim > cutoff else False

    @staticmethod
    def is_alignment(align):

        align = AlignIO.read(align, "fasta")

        return True if align else False
