import pyviennacl as p
import numpy as np


def _solve_cuda(lap_sparse, B, return_full_prob=False, maxiter=100, tol=5e-5):
    # get row col data
    lap_size = lap_sparse.get_shape()
    lap_sparse = lap_sparse.tocoo()
    rows = lap_sparse.row
    cols = lap_sparse.col
    data = lap_sparse.data
    lap_p = p.CompressedMatrix(lap_size[0], lap_size[1])
    for row, col, val in zip(rows, cols, data):
        row = int(row)
        col = int(col)
        lap_p[row, col] = float(val)
    X = []
    for i in range(len(B)):
        bi = -B[i].todense()
        bi = np.asarray(bi)
        bi = bi.reshape(-1)
        bi = bi.astype('float64')
        B_p = p.Vector(bi)
        x = p.solve(lap_p, B_p, p.upper_tag())
        X.append(x)
    if not return_full_prob:
        X = np.array(X)
        X = np.argmax(X, axis=0)
    return X