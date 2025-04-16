"""
@Env: /anaconda3/python3.11
@Time: 2025/4/16-10:04
@Auth: karlieswift
@File: eev.py
@Desc: 
"""

import numpy as np
import itertools


class EnergyEntropy():
    def __init__(self, data_type="DNA", energy_values=2, mutual_information_energy=2):
        """
        data_type (str):  The data type can be "protein", "DNA" or "RNA"
        energy_values (int):  The incoming energy value: 1,2,3,4.....
        mutual_information_energy (int):  The energy of mutual information : 2,3,4.....
        """
        super(EnergyEntropy, self).__init__()
        self.data_type = data_type
        self.energy_values = energy_values
        self.mutual_information_energy = mutual_information_energy
        self.data_type_elements = self._get_data_type(data_type=data_type)

    def seq2vector(self, seq):
        number_XY = self._get_pairs(self.data_type_elements)
        seq_len = len(seq)
        for i in range(seq_len):
            s1 = seq[i:i + 2]
            if s1 in number_XY:
                number_XY[s1] += 1

        # E1
        number_X = {e: seq.count(e) for e in self.data_type_elements}
        # 统计AGCT的概率
        p_X = {e: number_X[e] / seq_len for e in self.data_type_elements}

        H_X = {k: 0 for k, v in p_X.items()}
        for k, v in p_X.items():
            if v != 0:
                H_X[k] = -v * np.log2(v)

        E1 = []

        for energy_n in range(1, self.energy_values):
            E1 += self._compute_combinations_expand_info(H_X, number_X, k=energy_n)

        H_p_point_AGCT = {}
        for e in self.data_type_elements:
            temp = {k: v for k, v in number_XY.items() if k[1] == e}
            temp = {k: v / sum(number_XY.values()) for k, v in temp.items()}  # energyEntropymax
            H_p_point_AGCT[e] = sum([-v * np.log2(v + 1e-10) for v in temp.values()])

        E2 = []
        for energy_n in range(1, self.energy_values):
            E2 += self._compute_combinations_expand_info(H_p_point_AGCT, number_X, k=energy_n)

        # E3:
        all_postion = (seq_len + 1) * seq_len * 0.5

        postion_number = {e: [] for e in self.data_type_elements}
        for i in range(seq_len):
            if seq[i] in self.data_type_elements:
                postion_number[seq[i]].append(i + 1)

        postion_number_merge = {k: sum(v) for k, v in postion_number.items()}
        relative_postion = {k: all_postion - v for k, v in postion_number_merge.items()}
        relative_postion_p = {k: v / all_postion for k, v in relative_postion.items()}
        H_relative_postion_p = {k: 0 for k in self.data_type_elements}
        for k, v in relative_postion_p.items():
            if v != 0:
                H_relative_postion_p[k] = -v * np.log2(v)

        # E3
        E3 = []
        for energy_n in range(1, self.energy_values):
            E3 += self._compute_combinations_expand_info(H_relative_postion_p, number_X, k=energy_n)

        # E4
        E4 = self._get_mutual_information_multiple(seq, r=self.mutual_information_energy)
        all_vector = E1 + E2 + E3 + E4

        return [float(v) for v in all_vector]

    def _get_data_type(self, data_type):
        data_type_dict = {
            "dna": 'ACGT',
            "protein": 'ACDEFGHIKLMNPQRSTVWY',
            "RNA": 'ACGU',
        }
        return data_type_dict.get(data_type.lower())

    def _compute_combinations_expand_info(self, I_XY, number_XY, k):
        keys = I_XY.keys()
        combinations = itertools.combinations(keys, k)
        results = {}
        for comb in combinations:
            # 计算 I_XY 的乘积
            product_I_XY = 0
            product_number_XY = 0
            for key in comb:
                product_I_XY += I_XY[key]
                product_number_XY += number_XY[key]
            results[''.join(comb)] = product_I_XY * product_number_XY
        return list(results.values())

    def _get_mutual_information_multiple(self, seq, r):
        mutual_information_number = self._generate_combinations(elements=self.data_type_elements, r=r)
        # print(mutual_information_number)
        mutual_information_number_dict = {k: 0 for k in mutual_information_number}
        p_i = {e: seq.count(e) / len(seq) for e in self.data_type_elements}
        p_x_y_z = {k: 0 for k in mutual_information_number}
        for x_y_z in mutual_information_number:
            p_x_y_z[x_y_z] = p_i[x_y_z[0]]
            for x_y_z_i in x_y_z[1:]:
                p_x_y_z[x_y_z] *= p_i[x_y_z_i]

        for i in range(len(seq) - (r - 1)):
            s = "".join(sorted(seq[i:i + r]))
            if s in mutual_information_number:
                mutual_information_number_dict[s] = mutual_information_number_dict.get(s, 0) + 1
        P_mutual_information_number_dict = {k: v / len(seq) for k, v in mutual_information_number_dict.items()}
        # energy_info={k:mutual_information_number_dict[k]*v*np.log2(v/p_x_y_z[k]) for k,v in P_mutual_information_number_dict}
        energy_info = {}
        for k, v in P_mutual_information_number_dict.items():
            if v == 0:
                energy_info[k] = 0
            else:
                energy_info[k] = np.log2(v / p_x_y_z[k]) * v * mutual_information_number_dict[k]

        return list(energy_info.values())

    def _get_pairs(self, elements):
        combinations = list(itertools.product(elements, repeat=2))
        combinations = {a + b: 0 for (a, b) in combinations}

        return combinations

    def _generate_combinations(self, elements, r):
        return ["".join(sorted(comb)) for comb in itertools.combinations(elements, r)]


if __name__ == '__main__':
    # a=entropy_encoder2('EVQVVESGGGAVQPGGSLTLSCAASEPVFEANTMGWYRQAPGKQRELVASIRNVGGTNYADSVKGRFTITRGAAKNTINLQMNDLKPEDTAVYYCNAVIQLLFDSREYWGQGTQVTVSS')
    # seq="AAACA"
    seq = "GGTCAATTTAAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTATCGAATCCGA"
    # seq="EVQVVESGGGAVQPGGSLTLSCAASEPVFEANTMGWYRQAPGKQRELVASIRNVGGTNYADSVKGRFTITRGAAKNTINLQMNDLKPEDTAVYYCNAVIQLLFDSREYWGQGTQVTVSS"
    seq = "AAGATTCGGGCCTTCGGGTCCACCCGTTCCGCAGCTGTGCGCTCTTTGGGCTGCACGCTGTGTGATACACAACCCTCACACCTGTGAACGTATCGGGGGCGCGTAAGCGCTTCTGCTCAAAACATTTAACTACTTATGTTCAGAATGTAAAAAACTATAACAAATAACAACTTTCAACAACGGATCTCTTGGCTCTCGCATCGATGAAGAACGCAGCGAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACCTTGCGCTCTGTGGTATTCCGCAGAGCATGCCTGTTTGAGTGTCACGTAAACCATCGCCCTTGGGATTTCGATCTCTATGAGGTGGACTTGGACTGTGCCGTAGTCGGCTCGTCTTGAAATGAATTAGCTTGCGCTCTTTAGAGTGTCCGGCACCGGTGTGATAATTATCTGCGCCAACGCCTATGGCCTCTTCTTGCGGTGCTGCTTACAGTAGTCCGAAAGGACAGATCTACTTTAAAGCTTTGGCCTCA"
    # seq="AAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCGAAAACAAAAGCG"
    seq = "GAGTAAGATC"
    eev = EnergyEntropy()
    vector = eev.seq2vector(seq)

    print(len(vector))
    print(vector)
