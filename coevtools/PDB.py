from Bio.PDB import Polypeptide, PDBParser, is_aa, Selection, NeighborSearch
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceCalculator


class PDB:

    """" Parse PDB file """

    def __init__(self, pdb_file, ref="cb"):
        """
        Create PDB object

        Parameters:
        -------------------
        pdb_file : Atomic coordinates from protein complex.
        ref: Reference for amino acids where 'cb' is beta carbon atom and 'center' is geometric center.

        """
        if ref != "cb" and ref != "center":
            raise ValueError("Reference must be 'cb' or 'center'.")

        self.__pdb_file = pdb_file
        self.__ref = ref
        self.__name = ""
        self.__cutoff = 8
        self.__res_list = [[], []]
        self.__res_id_list = []
        self.__pairs_by_cb = []
        self.__pairs_by_center = []
        self.__parser()
        self.__get_ids_list()

    def __reference(self):

        if self.__ref == "cb":
            if not self.__pairs_by_cb:
                self.__contacts_by_cb()
            return self.__pairs_by_cb
        elif self.__ref == "center":
            if not self.__pairs_by_center:
                self.__contacts_by_center()
            return self.__pairs_by_center

    def __centroid(self, res):

        coords = np.zeros(3)
        atoms = Selection.unfold_entities(res, 'A')
        for atom in atoms:
            coords += atom.get_coord()

        return coords / len(atoms)

    def __parser(self):

        pdb = PDBParser(QUIET=True)
        structure = pdb.get_structure(self.__pdb_file, self.__pdb_file)
        self.__model = structure[0]

        chain_res_list = [Selection.unfold_entities(chain, "R") for chain in self.__model]
        for i, chain in enumerate(chain_res_list):
            for aa in chain:
                if is_aa(aa.get_resname(), standard=True):
                    self.__res_list[0].append(aa) if i == 0 else self.__res_list[1].append(aa)

    def __get_ids_list(self):


        for chain in self.__res_list:
            res_id = []
            for residue in chain:
                if is_aa(residue.get_resname(), standard=True):
                    res_id.append(residue.get_id()[1])
            self.__res_id_list.append(res_id)

    def __ref_atoms(self):

        atom_list = []

        for chain_res_list in self.__res_list:
            for residue in chain_res_list:
                if is_aa(residue):
                    if residue.has_id("CB"):
                        atom_list.append(residue["CB"])
                    else:
                        atom_list.append(residue["CA"])

        return atom_list

    def __missing_ids(self):

        missing_ids = []

        for idx in self.__res_id_list:
            start, end = idx[0], idx[-1]
            missing_ids.append(sorted(set(range(start, end + 1)).difference(idx)))

        return missing_ids

    def __contacts_by_cb(self):

        atom_list = self.__ref_atoms()

        ns = NeighborSearch(atom_list)

        for atom in atom_list:
            if atom.get_parent() in self.__res_list[0]:
                neighbors = ns.search(atom.get_coord(), self.__cutoff, "R")

                for res in neighbors:
                    if res not in self.__res_list[0]:
                        self.__pairs_by_cb.append((atom.get_parent(), res))

    def __contacts_by_center(self):

        for res_a in self.__res_list[0]:
            res_a_coord = self.__centroid(res_a)
            for res_b in self.__res_list[1]:
                res_b_coord = self.__centroid(res_b)
                if np.sqrt(np.sum(np.float_power(res_a_coord - res_b_coord, 2))) <= self.__cutoff:
                    self.__pairs_by_center.append([res_a, res_b])

    @property
    def fasta(self):

        fasta_res_list = []
        for chain in self.__res_list:
            fasta_chain = []
            for res in chain:
                if is_aa(res.get_resname(), standard=True):
                    fasta_chain.append(Polypeptide.three_to_one(res.get_resname()))

            fasta_chain = ''.join(fasta_chain)
            fasta_res_list.append(fasta_chain)

        return fasta_res_list

    #TODO
    def correct_shift(self, residue, add_a_ids=False):

        count = len(self.__res_list[0]) if add_a_ids else 0

        a_shift = self.__res_list[0][0].get_id()[1]
        b_shift = self.__res_list[1][0].get_id()[1]

        res_id = 0
        missing_ids = self.__missing_ids()

        if residue in self.__res_list[0]:
            if missing_ids[0]:
                for i in missing_ids[0]:
                    if residue.get_id()[1] > i:
                        res_id -= 1
            res_id += residue.get_id()[1] - a_shift

        else:
            if missing_ids[1]:
                for i in missing_ids[1]:
                    if residue.get_id()[1] > i:
                        res_id -= 1
            res_id += residue.get_id()[1] - b_shift + count

        return res_id

    @property
    def ICs(self):

        pairs = self.__reference()
        flat_list = [res for pair in pairs for res in pair]
        ics = sorted(list(set([self.correct_shift(ic, add_a_ids=True) for ic in flat_list])))
        return ics

    @property
    def pairs_shift(self):

        pairs = self.__reference()
        c_map = np.asarray([[self.correct_shift(id, add_a_ids=True) for id in pair] for pair in pairs])
        return c_map

    @property
    def pairs(self):

        pairs = self.__reference()
        c_map = np.asarray([list(map(self.correct_shift, pair)) for pair in pairs])
        return c_map

    @property
    def A_ICs(self):

        pairs = self.__reference()
        a_ics = [i[0] for i in pairs]
        a_ics = sorted(list(set(map(self.correct_shift, a_ics))))
        return a_ics

    @property
    def shift_B_ICs(self):

        pairs = self.__reference()
        b_ics = [i[1] for i in pairs]
        b_ics = sorted(list(set([self.correct_shift(ic, add_a_ids=True) for ic in b_ics])))
        return b_ics

    @property
    def B_ICs(self):

        pairs = self.__reference()
        b_ics = [i[1] for i in pairs]
        b_ics = sorted(list(set(map(self.correct_shift, b_ics))))
        return b_ics

    @property
    def pairs_wtt_correction(self):

        pairs = self.__reference()
        c_map = [list(map(lambda x: x.get_id()[1], pair)) for pair in pairs]
        return c_map

    @property
    def all_ICs(self):

        return [list(range(len(self.__res_list[0]))),
                [x + len(self.__res_list[0]) for x in list(range(len(self.__res_list[1])))]]

    @property
    def n_pairs(self):

        pairs = self.__reference()
        return len(pairs)

    @property
    def id_list(self):
        return self.__res_id_list

    @property
    def ref(self):
        return self.__ref

    @ref.setter
    def ref(self, value):
        self.__ref = value

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__name = name

    @property
    def cutoff(self):
        return self.__cutoff

    @cutoff.setter
    def cutoff(self, value):
        self.__cutoff = value

