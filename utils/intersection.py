"""
@Env: /anaconda3/python3.11
@Time: 2025/2/7-16:26
@Auth: karlieswift
@File: intersection.py
@Desc: 
"""

import numpy as np
from scipy.optimize import linprog


def intersection(mutset1, mutset2):
    m, l = mutset1.shape
    n, _ = mutset2.shape

    c = np.ones((m + n,))

    A0 = np.hstack((mutset1.T, -mutset2.T))

    a1 = np.hstack((np.ones((1, m)), np.zeros((1, n))))
    b1 = np.hstack((np.zeros((1, m)), np.ones((1, n))))
    Aeq = np.vstack((A0, a1, b1))

    beq = np.hstack((np.zeros(l), 1, 1))

    lb = np.zeros((m + n,))
    ub = np.ones((m + n,))


    res = linprog(c, A_ub=None, bounds=list(zip(lb, ub)), A_eq=Aeq, b_eq=beq, method='highs')

    if res.success:
        T = 0
    else:
        T = 1

    return T
if __name__ == '__main__':
    mutset1 = np.random.rand(100, 10)
    mutset2 = np.random.rand(100, 10) + 0.5
    result = intersection(mutset1, mutset2)
    print(result)


