"""
@Env: /anaconda3/python3.11
@Time: 2025/4/16-10:34
@Auth: karlieswift
@File: cnv.py
@Desc: 
"""


"""
@Env: /anaconda3/python3.11
@Time: 2025/4/16-10:33
@Auth: karlieswift
@File: nv.py
@Desc: 
"""


import numpy as np
import itertools

elements = 'ACGT'


# elements ='ACDEFGHIKLMNPQRSTVWY'

def get_pairs():
    combinations = list(itertools.product(elements, repeat=2))
    combinations = {a + b: 0 for (a, b) in combinations}
    return combinations



def calculate_covariance(sequence, mus, ns):
    unique_nucleotides = sorted(set(sequence))
    pairs = [''.join(pair) for pair in itertools.combinations(unique_nucleotides, 2)]
    len_seq = len(sequence)
    positions = np.arange(1, len_seq + 1)
    indicators = {pair: np.array([1 if x in pair else 0 for x in sequence]) for pair in pairs}
    covariances = {}
    cov_val = []
    for pair, indicator in indicators.items():
        k, l = pair[0], pair[1]
        if k not in elements or l not in elements:
            continue
        cov = np.sum((positions - mus[k]) * (positions - mus[l]) * indicator)
        cov /= (ns[k] * ns[l])
        # cov /= (np.sqrt(ns[k]) * np.sqrt(ns[l]) * len_seq)
        cov_val.append(cov)
        covariances[pair] = cov
    # print(covariances)

    # all_6pair_cov = {'AC': 0, 'AG': 0, 'AT': 0, 'CG': 0, 'CT': 0, 'GT': 0}
    all_6pair_cov = {a[0] + a[1]: 0 for a in list(itertools.combinations(elements, 2))}
    # print(len(all_6pair_cov))
    # Update the all_6pair_cov dictionary with the given covariances
    for pair in covariances:
        if pair in all_6pair_cov:
            all_6pair_cov[pair] = covariances[pair]

    return list(all_6pair_cov.values())


def calculate_values(sequence):
    n = {nucl: 0 for nucl in elements}
    mu = {nucl: 0 for nucl in elements}
    D2 = {nucl: 0 for nucl in elements}

    for nucl in elements:
        n[nucl] = sequence.count(nucl)

    for nucl in elements:
        positions = [i + 1 for i, char in enumerate(sequence) if char == nucl]
        if positions:
            mu[nucl] = sum(positions) / n[nucl]

    total_length = len(sequence)
    for nucl in elements:
        positions = [i + 1 for i, char in enumerate(sequence) if char == nucl]
        if positions:
            D2[nucl] = sum((pos - mu[nucl]) ** 2 for pos in positions) / total_length / n[nucl]

    return n, mu, D2


def calculate_NV_12(sequence):
    n, mu, D2 = calculate_values(sequence)
    NV_8 = list(n.values()) + list(mu.values())  # NV_8 = [nA, nG, nC, nT, muA, muG, muC, muT]
    NV_12 = NV_8 + list(D2.values())  # D2_4 = [D2A, D2G, D2C, D2T]
    return NV_12


def calculate_NV_18(sequence):
    n, mu, D2 = calculate_values(sequence)
    NV_8 = list(n.values()) + list(mu.values())  # NV_8 = [nA, nG, nC, nT, muA, muG, muC, muT]
    NV_12 = NV_8 + list(D2.values())  # D2_4 = [D2A, D2G, D2C, D2T]
    NV_6 = calculate_covariance(sequence, mu, n)
    NV_18 = NV_12 + NV_6
    return [float( i) for i in NV_18]


if __name__ == '__main__':
    seq = 'AGCTAAGA'

    print(len(calculate_NV_18(seq)))
    # print(len(entropy_markov_all(seq)))
    # print(joint_entropy(seq))
