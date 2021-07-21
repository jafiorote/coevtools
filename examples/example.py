from coevtools.PDB import PDB
from coevtools.MSA import MSA
from coevtools.Information import Information
import numpy as np

# parse pdb file
pdb = PDB("1BXR_AB.pdb")

# parse msa fasta files
msa = MSA("1BXR_A.fas", "1BXR_B.fas")

# Shannon information
site_freq = Information.site_freq(msa, pdb.ICs)
pair_freq = Information.pair_freq(msa, pdb.ICs, site_freq)
mi_arr, h_arr = Information.shannon_info(msa, pdb, site_freq, pair_freq)

print(f"Total mutual information = {np.sum(mi_arr)}")
print(f"Total entropy = {np.sum(h_arr)}")