from __future__ import division
import pycuda.autoinit
import pycuda.driver as drv
import pycuda.gpuarray as gpuarray
from pycuda.tools import DeviceMemoryPool, PageLockedMemoryPool
from pycuda.sparse.packeted import PacketedSpMV
from pycuda.sparse.operator import DiagonalPreconditioner
from pycuda.sparse.cg import solve_pkt_with_cg

import numpy as np

def _solve_cuda(lap_sparse, B, return_full_prob=False, maxiter=100, tol=5e-5):
    """
    solves lap_sparse X_i = B_i for each phase i, using the conjugate
    gradient method. For each pixel, the label i corresponding to the
    maximal X_i is returned.
    """
    print("using gpu mode")
    dev_pool = DeviceMemoryPool()
    pagelocked_pool = PageLockedMemoryPool()
    csr_mat = lap_sparse
    csr_mat = csr_mat.astype(np.float32)
    inv_mat_diag = 1 / csr_mat.diagonal()
    spmv = PacketedSpMV(csr_mat, True, csr_mat.dtype)
    X = []
    for i in range(len(B)):
        rhs = -B[i].astype(spmv.dtype)
        if True:
            precon = DiagonalPreconditioner(spmv.permute(gpuarray.to_gpu(inv_mat_diag,
                                                                         allocator=dev_pool.allocate)))
        else:
            precon = None
        print("start solve")
        start = drv.Event()
        stop = drv.Event()
        start.record()
        rhs_gpu = gpuarray.to_gpu(rhs, dev_pool.allocate)
        tol = 1e-7 if spmv.dtype == np.float64 else tol
        res_gpu, it_count, res_count = solve_pkt_with_cg(spmv, rhs_gpu,
                                                         precon, tol=tol,
                                                         pagelocked_allocator=pagelocked_pool.allocate)
        res = res_gpu.get()
        stop.record()
        stop.synchronize()
        elapsed = stop.time_since(start) * 1e-3
        est_flops = (csr_mat.nnz * 2 * (it_count + res_count)
                     + csr_mat.shape[0] * (2 + 2 + 2 + 2 + 2) * it_count)
        if precon is not None:
            est_flops += csr_mat.shape[0] * it_count
        print("size: %d, elapsed: %g s, %d it, %d residual, it/second: %g, "
              "%g gflops/s" % (
                  csr_mat.shape[0],
                  elapsed, it_count, res_count, it_count / elapsed,
                  est_flops / elapsed / 1e9))
        x0 = res[0]
        X.append(x0)
    pagelocked_pool.stop_holding()
    dev_pool.stop_holding()
    if not return_full_prob:
        X = np.array(X)
        X = np.argmax(X, axis=0)
    return X
