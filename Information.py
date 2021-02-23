import numpy as np
from scipy import spatial
from Bio.Phylo.TreeConstruction import DistanceCalculator

class Information:

    @classmethod
    def sequence_weight(cls, msa, theta):

        hamming_distance = spatial.distance.pdist(msa.full_encoded, "hamming")
        weight_matrix = spatial.distance.squareform(hamming_distance < (1.0 - theta))
        return 1.0 / (np.sum(weight_matrix, axis=1) + 1.0)

    @classmethod
    def site_freq(cls, msa, ics, lambda_=0.001):

        site_freq = np.zeros((len(ics), msa.q), dtype=float)

        for i in range(len(ics)):
            vec = np.bincount(msa.sub_full_MSA(ics)[:, i], weights=msa.sequence_weight)
            site_freq[i, 0:vec.size] = vec

        site_freq /= msa.Meff
        return (1 - lambda_) * site_freq + lambda_ / msa.q

    @classmethod
    def pair_freq(cls, msa, ics, sitefreq_array, lambda_=0.001, swap_idx=None):

        pair_freq = np.zeros((len(ics), msa.q, len(ics), msa.q), dtype=float)
        encoded = msa.sub_full_MSA(ics, swap_idx)
        for i in range(len(ics)):
            for j in range(len(ics)):
                np.add.at(pair_freq, [i, encoded[:, i], j, encoded[:, j]], msa.sequence_weight.reshape(-1, 1)[:, 0])

        pair_freq /= msa.Meff
        pair_freq = (1 - lambda_) * pair_freq + lambda_ / (msa.q * msa.q)

        for i, aa in enumerate(ics):
            for am_i in range(msa.q):
                for am_j in range(msa.q):
                    if am_i == am_j:
                        pair_freq[i, am_i, i, am_j] = sitefreq_array[i, am_i]
                    else:
                        pair_freq[i, am_i, i, am_j] = 0.0

        return pair_freq

    @classmethod
    def shannon_info(cls, msa, model, site_freq, pair_freq):

        contact_map = model.contact_map
        A_ICs = model.A_ICs
        B_ICs = model.shift_B_ICs
        nA = len(A_ICs)
        nB = len(B_ICs)

        mi_matrix_pp = np.zeros((nA, nB), dtype=float)
        h_matrix = np.zeros((nA, nB), dtype=float)

        for i, col_i in enumerate(A_ICs):
            for j, col_j in enumerate(B_ICs):
                tmp1 = (np.sum(pair_freq[i, :, j + nA, :] *
                               np.log(pair_freq[i, :, j + nA, :] /
                               (np.transpose(np.broadcast_to(site_freq[i, :], (msa.q, msa.q))) *
                                np.broadcast_to(site_freq[j + nA, :], (msa.q, msa.q))))))
                tmp2 = (np.sum(pair_freq[i, :, j + nA, :] * np.log(pair_freq[i, :, j + nA, :]))) * -1
                if any((x == [col_i, col_j]).all() for x in contact_map):
                    mi_matrix_pp[i, j] = tmp1
                    h_matrix[i, j] = tmp2

        return mi_matrix_pp, h_matrix

    @classmethod
    def tsallis(cls, sitefreq, pairfreq, pairs, idx_pairs, q=21, alpha=1.75):

        nP = len(pairs)
        h_matrix_pp = np.zeros((nP, nP), dtype=float)
        mi_matrix_pp = np.zeros((nP, nP), dtype=float)

        for k in range(len(pairs)):
            i = idx_pairs[k][0]
            j = idx_pairs[k][1]
            pij = pairfreq[i, :, j, :]
            pi = np.transpose(np.broadcast_to(sitefreq[i, :], (q, q)))
            pj = np.broadcast_to(sitefreq[j, :], (q, q))
            h_matrix_pp[i, j] = -np.sum((pij ** alpha) * ((pij ** (1 - alpha) - 1) / (1 - alpha)))
            mi_matrix_pp[i, j] = np.sum((pij ** alpha) * (((pij / (pi * pj)) ** (1 - alpha) - 1) / (1 - alpha)))

        return mi_matrix_pp, h_matrix_pp

    @classmethod
    def mirror_tree(cls, msaA, msaB):

        calculator = DistanceCalculator('blosum62')
        dm_A = calculator.get_distance(msaA)
        dm_B = calculator.get_distance(msaB)
        av_A = np.average(dm_A)
        av_B = np.average(dm_B)

        num = 0
        den_A = 0
        den_B = 0
        num_of_seqs = len(msaA)

        for i in range(num_of_seqs):
            for j in range(num_of_seqs):
                num += (dm_A[i, j] - av_A) * (dm_B[i, j] - av_B)
                den_A += (dm_A[i, j] - av_A) ** 2
                den_B += (dm_B[i, j] - av_B) ** 2
        r = num / (np.sqrt(den_A) * np.sqrt(den_B))

        # corr = np.corrcoef(dm_A, den_B)

        return r

