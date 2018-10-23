#encoding: utf-8
"""
Time：
    2018-4-25
Brief：
    python code for confidence estimation using random walks.
    support 2d image and 3d volume
Ref1:
    Ultrasound confidence maps
    A. Karamalis, W. Wein, T. Klein, N. Navab: Ultrasound Confidence Maps
    using Random Walks, Medical Image Analysis, 16, 6, 1101 - 1112, 2012
    DOI: http://dx.doi.org/10.1016/j.media.2012.07.005
Ref2:
    skimage\segmentation\random_walker_segmentation.py
"""
import numpy as np
from scipy.signal import hilbert2, hilbert
from scipy import sparse, ndimage as ndi
from scipy.sparse.linalg import cg
from skimage import img_as_float
# from skimage.filters import rank_order
from skimage.filters import rank_order
from cuda_solver import _solve_cuda
# from pycuda_solver import _solve_cuda

try:
    from scipy.sparse.linalg.dsolve import umfpack
    old_del = umfpack.UmfpackContext.__del__

    def new_del(self):
        try:
            old_del(self)
        except AttributeError:
            pass
    umfpack.UmfpackContext.__del__ = new_del
    UmfpackContext = umfpack.UmfpackContext()
except:
    UmfpackContext = None

try:
    from pyamg import ruge_stuben_solver
    amg_loaded = True
except ImportError:
    amg_loaded = False


def normalize_data(data, eps=1e-7):
    data = (data - np.min(data[:])) / ((np.max(data[:]) - np.min(data[:])));
    return data


def convert_data_to_4d(data, multichannel):
    # convert data to 4d
    if not multichannel:
        if data.ndim < 2 or data.ndim > 3:
            raise ValueError('single-channel image and volume must have 2 or 3 dims')
        dims = data.shape  # To reshape final labeled result
        # convert 2d image to 4d
        data = np.atleast_3d(img_as_float(data))[..., np.newaxis]
    else:
        if data.ndim < 3:
            raise ValueError('multi-channel image and volume must have 3 or 4 dims')
        dims = data[..., 0].shape  # To reshape final labeled result
        data = img_as_float(data)
        if data.ndim == 3:          # 2D multispectral, needs singleton in 3rd axis
            data = data[:, :, np.newaxis, :]
    return data, dims

def spacing_check(spacing, dims):
    # Spacing kwarg checks
    if spacing is None:
        spacing = np.asarray((1.,) * 3)
    elif len(spacing) == len(dims):
        if len(spacing) == 2:  # Need a dummy spacing for singleton 3rd dim
            spacing = np.r_[spacing, 1.]
        else:                  # Convert to array
            spacing = np.asarray(spacing)
    else:
        raise ValueError('Input argument `spacing` incorrect, should be an '
                         'iterable with one number per spatial dimension.')
    return spacing

def labels_check(labels):
    labels = np.copy(labels)
    label_values = np.unique(labels)
    # Reorder label values to have consecutive integers (no gaps)
    if np.any(np.diff(label_values) != 1):
        mask = labels >= 0
        labels[mask] = rank_order(labels[mask])[0].astype(labels.dtype)
    labels = labels.astype(np.int32)
    if np.any(labels < 0):
        filled = ndi.binary_propagation(labels > 0, mask=labels >= 0)
        labels[np.logical_and(np.logical_not(filled), labels == 0)] = -1
        del filled
    labels = np.atleast_3d(labels)
    return labels


def attenuation_weighting(data_size, alpha):
    """
    Compute attenuation weighting, defalt attenuate along first axis
    Args:
        data_size: sequence, size of ultrasound image, (height, width, depth)
                   do not include channel axis
        alpha: float, Attenuation coefficient according to Beer–Lambert law
    Return:
        W: Weighting expresing depth-dependent attenuation
    """
    depth = data_size[0]
    dw = np.linspace(1, depth, depth)
    dw = dw / depth
    width = data_size[1]
    if len(data_size) == 2:
        dw = np.reshape(dw, (depth, 1, 1, 1))
        dw = np.tile(dw, [1, width, 1, 1])
    elif len(data_size) == 3:
        dw = np.reshape(dw, (depth, 1, 1, 1))
        frames = data_size[2]
        dw = np.tile(dw, [1, width, frames, 1])
    dw = normalize_data(dw)
    dw = 1.0 - np.exp(-alpha * dw)

    return dw


def _clean_labels_ar(X, labels, copy=False):
    X = X.astype(labels.dtype)
    if copy:
        labels = np.copy(labels)
    labels = np.ravel(labels)
    labels[labels == 0] = X
    return labels


def _make_graph_edges_3d(n_x, n_y, n_z):
    """Returns a list of edges for a 3D image.
    Args:
        n_x: integer
            The size of the grid in the x direction.
        n_y: integer
            The size of the grid in the y direction
        n_z: integer
            The size of the grid in the z direction
    Returns:
        edges : (2, N) ndarray
            with the total number of edges::
                N = n_x * n_y * (nz - 1) +
                    n_x * (n_y - 1) * nz +
                    (n_x - 1) * n_y * nz
            Graph edges with each column describing a node-id pair.
    """
    vertices = np.arange(n_x * n_y * n_z).reshape((n_x, n_y, n_z))
    edges_deep = np.vstack((vertices[:, :, :-1].ravel(),
                            vertices[:, :, 1:].ravel()))
    edges_right = np.vstack((vertices[:, :-1].ravel(),
                             vertices[:, 1:].ravel()))
    edges_down = np.vstack((vertices[:-1].ravel(), vertices[1:].ravel()))
    edges = np.hstack((edges_deep, edges_right, edges_down))
    return edges


def _compute_gradients_3d(data, spacing):
    gr_deep = np.abs(data[:, :, :-1] - data[:, :, 1:]).ravel() / spacing[2]
    gr_right = np.abs(data[:, :-1] - data[:, 1:]).ravel() / spacing[1]
    gr_down = np.abs(data[:-1] - data[1:]).ravel() / spacing[0]

    return np.r_[gr_deep, gr_right, gr_down]


def _compute_weights_3d(data, spacing, gamma, beta=130, eps=1.e-6, multichannel=False):
    """s
    creat weighted matrix
    Weight calculation is main difference in multispectral version
    Original gradient**2 replaced with sum of gradients ** 2
    Args:
        gamma: float, horizontal penalty
    """
    gradients = 0
    for channel in range(0, data.shape[-1]):
        # add horizental punish
        gradients += gamma
        # gradients += _compute_gradients_3d(data[..., channel],
        #                                    spacing) ** 2
        gradients += _compute_gradients_3d(data[..., channel],
                                           spacing)
    gradients = normalize_data(gradients)
    # All channels considered together in this standard deviation
    # beta /= 10 * data.std()
    # if multichannel:
    #     # New final term in beta to give == results in trivial case where
    #     # multiple identical spectra are passed.
    #     beta /= np.sqrt(data.shape[-1])
    print("beta value is ", beta)
    gradients *= beta
    weights = np.exp(- gradients)
    weights += eps
    return weights


def _mask_edges_weights(edges, weights, mask):
    """
    Remove edges of the graph connected to masked nodes, as well as
    corresponding weights of the edges.
    """
    mask0 = np.hstack((mask[:, :, :-1].ravel(), mask[:, :-1].ravel(),
                       mask[:-1].ravel()))
    mask1 = np.hstack((mask[:, :, 1:].ravel(), mask[:, 1:].ravel(),
                       mask[1:].ravel()))
    ind_mask = np.logical_and(mask0, mask1)
    edges, weights = edges[:, ind_mask], weights[ind_mask]
    max_node_index = edges.max()
    # Reassign edges labels to 0, 1, ... edges_number - 1
    order = np.searchsorted(np.unique(edges.ravel()),
                            np.arange(max_node_index + 1))
    edges = order[edges.astype(np.int64)]
    return edges, weights


def _make_laplacian_sparse(edges, weights):
    """
    Sparse implementation
    """
    pixel_nb = edges.max() + 1
    diag = np.arange(pixel_nb)
    i_indices = np.hstack((edges[0], edges[1]))
    j_indices = np.hstack((edges[1], edges[0]))
    data = np.hstack((-weights, -weights))
    lap = sparse.coo_matrix((data, (i_indices, j_indices)),
                            shape=(pixel_nb, pixel_nb))
    connect = - np.ravel(lap.sum(axis=1))
    lap = sparse.coo_matrix(
        (np.hstack((data, connect)), (np.hstack((i_indices, diag)),
                                      np.hstack((j_indices, diag)))),
        shape=(pixel_nb, pixel_nb))
    return lap.tocsr()


def _build_laplacian(data, spacing, gamma, mask=None, beta=50,
                     multichannel=False):
    """
    creat lattice
    Args:
        data: 4d array, volume or image
        spacing: 1d array, data spacing along every axis
        beta: float, random walk parm
        mask: 4d array, labeled mask
        gamma: float, horizontal walks penalty
    """
    l_x, l_y, l_z = tuple(data.shape[i] for i in range(3))
    edges = _make_graph_edges_3d(l_x, l_y, l_z)
    weights = _compute_weights_3d(data, spacing, gamma, beta=beta, eps=1.e-7,
                                  multichannel=multichannel)
    if mask is not None:
        edges, weights = _mask_edges_weights(edges, weights, mask)
    lap = _make_laplacian_sparse(edges, weights)
    del edges, weights
    return lap


def _buildAB(lap_sparse, labels):
    """
    Build the matrix A and rhs B of the linear system to solve.
    A and B are two block of the laplacian of the image graph.
    """
    labels = labels[labels >= 0]
    indices = np.arange(labels.size)
    unlabeled_indices = indices[labels == 0]
    seeds_indices = indices[labels > 0]
    # The following two lines take most of the time in this function
    B = lap_sparse[unlabeled_indices][:, seeds_indices]
    lap_sparse = lap_sparse[unlabeled_indices][:, unlabeled_indices]
    nlabels = labels.max()
    rhs = []
    for lab in range(1, nlabels + 1):
        mask = (labels[seeds_indices] == lab)
        fs = sparse.csr_matrix(mask)
        fs = fs.transpose()
        rhs.append(B * fs)
    return lap_sparse, rhs


def confidence_map(data, labels,
                   alpha, beta, gamma,
                   return_full_prob,
                   multichannel=False, 
                   spacing=None,
                   mode='cg',
                   tol=1.e-3):
    """
    1. if multichannel=True,be careful channel last.
    2. if voxels have physics spacing, please set spacing list,
        Anisotropic data is commonly encountered in medical imaging
    Args:
        data: numpy tensor, Image to be segmented in phases. Gray-level `data` can be two- or
            three-dimensional, [height, width, frames];
        multichannel: bool, multichannel data can be three- or four-
            dimensional (multichannel=True) with the highest dimension denoting channels.
        spacing:
            Data spacing is assumed isotropic unless the `spacing` vkeyword argument is used.
        labels: array of ints, of same shape as `data` without channels dimension 
            Array of seed markers labeled with different positive integers
            for different phases. Zero-labeled pixels are unlabeled pixels.
            Negative labels correspond to inactive pixels that are not taken
            into account (they are removed from the graph). If labels are not
            consecutive integers, the labels array will be transformed so that
            labels are consecutive. In the multiccgmghannel case, `labels` should have
            the same shape as a single channel of `data`, i.e. without the final
            dimension denoting channels. 
        tol: float, tolerance to achieve when solving the linear system
        return_full_prob : bool, default False,If True, the probability that a pixel belongs
            to each of the labels will be returned, instead of only the most likely label.
        spacing : list or tuple, Spacing between voxels in each spatial dimension. If `None`, then
            the spacing between pixels/voxels in each dimension is assumed 1.
    Return:
        output : ndarray
            * If `return_full_prob` is False, array of ints of same shape as
            `data`, in which each pixel has been labeled according to the marker
            that reached the pixel first by anisotropic diffusion.
            * If `return_full_prob` is True, array of floats of shape
            `(nlabels, data.shape)`. `output[label_nb, i, j]` is the probability
            that label `label_nb` reaches the pixel `(i, j)` first.
    """
    if (labels != 0).all():
        print('Random walker only segments unlabeled areas, where ',
             'labels == 0. No zero valued areas in labels were ',
             'found. Returning provided labels.')

        if return_full_prob:
            # Find and iterate over valid labels
            unique_labels = np.unique(labels)
            unique_labels = unique_labels[unique_labels > 0]

            out_labels = np.empty(labels.shape + (len(unique_labels),),
                                  dtype=np.bool)
            # set probility of labeled seeds to 1.0
            for n, i in enumerate(unique_labels):
                out_labels[..., n] = (labels == i)
        else:
            out_labels = labels
        return out_labels

    data, dims = convert_data_to_4d(data, multichannel)
    dw = attenuation_weighting(dims, alpha)
    data = data * dw
    spacing = spacing_check(spacing, dims)
    labels = labels_check(labels)
    if np.any(labels < 0):
        # 构建laplacian
        lap_sparse = _build_laplacian(data, spacing, gamma, mask=labels >= 0,
                                      beta=beta, multichannel=multichannel)
    else:
        lap_sparse = _build_laplacian(data, spacing, gamma, beta=beta, multichannel=multichannel)
    lap_sparse, B = _buildAB(lap_sparse, labels)
    # We solve the linear system
    # lap_sparse X = B
    # where X[i, j] is the probability that a marker of label i arrives
    # first at pixel j by anisotropic diffusion.
    if mode == 'cg':
        X = _solve_cg(lap_sparse, B, tol=tol,
                      return_full_prob=return_full_prob)
    if mode == 'cg_mg':
        if not amg_loaded:
            print("""pyamg (http://pyamg.org/)) is needed to use
                this mode, but is not installed. The 'cg' mode will be used
                instead.""")
            X = _solve_cg(lap_sparse, B, tol=tol,
                          return_full_prob=return_full_prob)
        else:
            X = _solve_cg_mg(lap_sparse, B, tol=tol,
                             return_full_prob=return_full_prob)
    if mode == 'bf':
        X = _solve_bf(lap_sparse, B,
                      return_full_prob=return_full_prob)
    if mode == 'gpu':
        X = _solve_cuda(lap_sparse, B, return_full_prob=return_full_prob)

    # Clean up results
    if return_full_prob:
        labels = labels.astype(np.float)
        X = np.array([_clean_labels_ar(Xline, labels, copy=True).reshape(dims)
                      for Xline in X])
        for i in range(1, int(labels.max()) + 1):
            mask_i = np.squeeze(labels == i)
            X[:, mask_i] = 0
            X[i - 1, mask_i] = 1
    else:
        X = _clean_labels_ar(X + 1, labels).reshape(dims)
    return X  


def _solve_cg(lap_sparse, B, tol, return_full_prob=False):
    """
    solves lap_sparse X_i = B_i for each phase i, using the conjugate
    gradient method. For each pixel, the label i corresponding to the
    maximal X_i is returned.
    """
    print("using cg mode tol---", tol)
    lap_sparse = lap_sparse.tocsc()
    X = []
    for i in range(len(B)):
        x0 = cg(lap_sparse, -B[i].todense(), tol=tol)[0]
        X.append(x0)
    if not return_full_prob:
        X = np.array(X)
        X = np.argmax(X, axis=0)
    return X


def _solve_cg_mg(lap_sparse, B, tol, return_full_prob=False, maxiter=100):
    """
    solves lap_sparse X_i = B_i for each phase i, using the conjugate
    gradient method with a multigrid preconditioner (ruge-stuben from
    pyamg). For each pixel, the label i corresponding to the maximal
    X_i is returned.
    """
    print('use cg mg mode, tol---', tol, 'max iterate---', maxiter)
    X = []
    ml = ruge_stuben_solver(lap_sparse)
    M = ml.aspreconditioner(cycle='V')
    for i in range(len(B)):
        x0 = cg(lap_sparse, -B[i].todense(), tol=tol, M=M, maxiter=maxiter)[0]
        X.append(x0)
    if not return_full_prob:
        X = np.array(X)
        X = np.argmax(X, axis=0)
    return X


def _solve_bf(lap_sparse, B, return_full_prob=False):
    """
    solves lap_sparse X_i = B_i for each phase i. An LU decomposition
    of lap_sparse is computed first. For each pixel, the label i
    corresponding to the maximal X_i is returned.
    """
    print("using bf mode")
    lap_sparse = lap_sparse.tocsc()
    solver = sparse.linalg.factorized(lap_sparse.astype(np.double))
    X = np.array([solver(np.array((-B[i]).todense()).ravel())
                  for i in range(len(B))])
    if not return_full_prob:
        X = np.argmax(X, axis=0)
    return X


def confidence_map2d(data, alpha=1.5, beta=90, gamma=0.03, spacing=None, data_mode='B', solver_mode="cg"):
    """2d confidence map

    # Args
        data: 3d numpy array, RF mode data，(height, width)
        mode: string, 'RF' or 'B' mode data
        alpha: float, distance(vertical) penalty
        beta: float, Random walks parameter
        gamma: float, horizontal penalty
    """
    data = data.astype('float')
    data = normalize_data(data)
    if data_mode == "RF":
        data = np.abs(hilbert2(data))
    labels = np.zeros_like(data)
    labels[0, :] = 1        # 探头元素
    labels[-1, :] = 2       # shadow元素
    conf_map = confidence_map(data, labels, alpha, beta, gamma, True, mode=solver_mode, spacing=spacing)
    conf_map = conf_map[0,:,:]
    return conf_map


def confidence_map3d(data, alpha=1.5, beta=90, gamma=0.03, spacing=None, data_mode='B', solver_mode="cg"):
    """3d confidence map

    # Args
        data: 3d numpy array, RF mode data，(height, width, depth)
        mode: string, 'RF' or 'B' mode data
        alpha: float, distance(vertical) penalty
        beta: float, Random walks parameter
        gamma: float, horizontal penalty
    """
    data = normalize_data(data)
    if data_mode == "RF":
        data = np.abs(hilbert(data))
    labels = np.zeros_like(data)
    labels[0, :, :] = 1             # elements of probe
    labels[-1:, :, :] = 2           # elements of shadow
    conf_map = confidence_map(data, labels, alpha, beta, gamma, True, spacing=spacing, mode=solver_mode)
    conf_map = conf_map[0,:, :, :]

    return conf_map


