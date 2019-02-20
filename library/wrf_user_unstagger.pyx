#!/usr/bin/python
import numpy as np
cimport numpy as np
cimport cython
from datetime import datetime

DTYPE = np.float32

def wrf_user_unstagger_ARW(varin, unstagDim):
    dims = ()
    dims = np.shape(varin)
    nd = np.shape(dims)[0]

    print(unstagDim, nd)
    if (unstagDim == "X") | (unstagDim == "U"):
        dimU = dims[nd - 1]
        if (nd == 5):
            varout = 0.5 * (
                varin[:, :, :, :, :dimU - 2] + varin[:, :, :, :, 1:dimU - 0])

        if (nd == 4):
            varout = 0.5 * (
                varin[:, :, :, :dimU - 1] + varin[:, :, :, 1:dimU - 0])

        if (nd == 3):
            varout = 0.5 * (varin[:, :, :dimU - 1] + varin[:, :, 1:dimU - 0])

        if (nd == 2):
            varout = 0.5 * (varin[:, :dimU - 1] + varin[:, 1:dimU - 0])

    if (unstagDim == "Y") | (unstagDim == "V"):
        dimV = dims[nd - 2]
        if (nd == 5):
            varout = 0.5 * (
                varin[:, :, :, :dimV - 1, :] + varin[:, :, :, 1:dimV - 0, :])

        if (nd == 4):
            varout = 0.5 * (
                varin[:, :, :dimV - 1, :] + varin[:, :, 1:dimV - 0, :])

        if (nd == 3):
            varout = 0.5 * (varin[:, :dimV - 1, :] + varin[:, 1:dimV - 0, :])

        if (nd == 2):
            varout = 0.5 * (varin[:dimV - 1, :] + varin[1:dimV - 0, :])

    if (unstagDim == "Z"):
        dimW = dims[nd - 3]
        if (nd == 5):
            varout = 0.5 * (
                varin[:, :, 0:dimW - 1, :, :] + varin[:, :, 1:dimW - 0, :, :])

        if (nd == 4):
            varout = 0.5 * (
                varin[:, 0:dimW - 1, :, :] + varin[:, 1:dimW - 0, :, :])

        if (nd == 3):
            varout = 0.5 * (varin[0:dimW - 1, :, :] + varin[1:dimW - 0, :, :])

        if (nd == 2):
            varout = 0.5 * (varin[:, 0:dimW - 1] + varin[:, 1:dimW - 0])


    if (unstagDim == "M"):
        print('staggared dim is "M", no support for this yet')
        varout = varin

    for ele in ["X", "U", "Y", "V", "Z", "M"]:
        if (unstagDim == ele):
            return varout

def unstagger_NMM_Z(varin):
    dims_in = np.shape(varin)
    n_dims = len(dims_in)
    dimW = dims_in[n_dims - 3]
    if (n_dims == 5):
        var_out = 0.5 * (
            varin[:, :, 0:dimW - 1, :, :] + varin[:, :, 1:dimW - 0, :, :])
    if (n_dims == 4):
        var_out = 0.5 * (
            varin[:, 0:dimW - 1, :, :] + varin[:, 1:dimW - 0, :, :])
    if (n_dims == 3):
        var_out = 0.5 * (varin[:, 0:dimW - 1] + varin[:, 1:dimW - 0])
    if (n_dims == 2):
        var_out = 0.5 * (varin[:, 0:dimW - 1] + varin[:, 1:dimW - 0])
    return var_out


@cython.boundscheck(False)
@cython.wraparound(False)
def get_gridmask_NMM(int dim_ns_out, int dim_we_out, str unstagDim):
    """
    This part is general and should not depend on the
    higher dimenions, essentially setting up how the 
    grid is structured
    """

    begin = datetime.now()
    cdef Py_ssize_t dim_ns = dim_ns_out
    cdef Py_ssize_t dim_we = dim_we_out
    grid_mask = np.zeros((dim_ns, dim_we), dtype=DTYPE)
    cdef float[:, :] grid_mask_view = grid_mask

    if unstagDim == 'H':
        grid_mask_view[:, :] = np.nan
    if unstagDim == 'V':
        grid_mask_view[:, :] = 1
    if unstagDim == 'Z':
        grid_mask_view[:, :] = np.nan

    cdef Py_ssize_t ii, jj
    for ii in range(dim_ns):
        for jj in range(dim_we):
            if np.mod(ii, 2) == 0 and np.mod(jj, 2) == 0:
                if unstagDim == 'H' or unstagDim == 'Z':
                    grid_mask_view[ii, jj] = 1
                elif unstagDim == 'V':
                    grid_mask_view[ii, jj] = np.nan

            if np.mod(ii, 2) == 1 and np.mod(jj, 2) == 1:
                if unstagDim == 'H' or unstagDim == 'Z':
                    grid_mask_view[ii, jj] = 1
                if unstagDim == 'V':
                    grid_mask_view[ii, jj] = np.nan
    print(f'total time for grid_mask was {(datetime.now() - begin).total_seconds()}')
    return grid_mask

#@cython.boundscheck(False)
#@cython.wraparound(False)
#def regrid_NMM_rank2(np.ndarray[float, ndim=2] varin, np.ndarray[float, ndim=2] grid_mask, unstagDim):
#    """
#    Function to regrid the NMM data when dealing with
#    2-dimenional arrays.
#    Actually this probably never occurs due to the fact that if
#    we have 2-d data one of them is probably time or level
#    """
#
#    n_dims = len(np.shape(varin))
#    dim_ns_in = np.shape(varin)[n_dims - 2]
#    dim_we_in = np.shape(varin)[n_dims - 1]
#    dim_ns_out = dim_ns_in
#    dim_we_out = dim_we_in * 2
#
#
#    # 2 times but no time, i dont think this ever occurs
#    var_out = np.zeros((dim_ns_out, dim_we_out), dtype=DTYPE)
#    var_out[:, :] = grid_mask[:, :]
#    ii_count = 0
#    jj_count = 0
#    for ii in range(dim_ns_out):
#        for jj in range(dim_we_out):
#            if not np.isnan(grid_mask[ii, jj]):
#                var_out[ii, jj] = varin[ii_count, jj_count]
#                jj_count = jj_count + 1
#                if jj_count == dim_we_in:
#                    jj_count = 0
#                    ii_count = ii_count + 1
#    # now use your copied values to fill in the blank spots
#    #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
#    ns_dim = np.shape(grid_mask)[0]
#    we_dim = np.shape(grid_mask)[1]
#    for ii in range(1,  ns_dim-1):
#        for jj in range(1, we_dim-1):
#            if np.isnan(grid_mask[ii, jj]):
#                var1 = var_out[ii - 1, jj]
#                var2 = var_out[ii + 1, jj]
#                var3 = var_out[ii, jj - 1]
#                var4 = var_out[ii, jj + 1]
#                var_out[ii, jj] = 0.25 * (var1 + var2 + var3 + var4)
#    # below fixes the "edge" problem caused by destaggering
#    return var_out[1:ns_dim-1, 1:we_dim-1]
#

@cython.boundscheck(False)
@cython.wraparound(False)
def regrid_NMM_rank3(np.ndarray[float, ndim=3] varin, np.ndarray[float, ndim=2] grid_mask, str unstagDim):
    """
    Function to regrid the NMM data when dealing with
    3-dimenional arrays
    """

    cdef float[:, :, :] varin_view = varin
    cdef int n_dims = len(np.shape(varin))
    cdef int dim_ns_in = np.shape(varin)[n_dims - 2]
    cdef int dim_we_in = np.shape(varin)[n_dims - 1]
    cdef int dim_ns_out = dim_ns_in
    cdef int dim_we_out = dim_we_in * 2
    cdef Py_ssize_t dim1 = np.shape(varin)[0]
    cdef Py_ssize_t dim_ns = dim_ns_out
    cdef Py_ssize_t dim_we = dim_we_out

    cdef float[:, :] grid_mask_view = grid_mask

    var_out = np.zeros((dim1, dim_ns, dim_we), dtype=DTYPE)
    cdef float[:, :, :] var_out_view = var_out
    var_out_view[:, :, :] = np.nan

    cdef Py_ssize_t ii, jj
    cdef Py_ssize_t ii_count = 0
    cdef Py_ssize_t jj_count = 0
    for ii in range(dim_ns):
        for jj in range(dim_we):
            if not np.isnan(grid_mask_view[ii, jj]):
                var_out_view[:, ii, jj] = varin_view[:, ii_count, jj_count]
                jj_count = jj_count + 1
                if jj_count == dim_we_in:
                    jj_count = 0
                    ii_count = ii_count + 1

    # now use your copied values to fill in the blank spots
    #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
    cdef Py_ssize_t ns_dim = np.shape(grid_mask)[0]
    cdef Py_ssize_t we_dim = np.shape(grid_mask)[1]
    cdef Py_ssize_t dim_1
    cdef float var1
    cdef float var2
    cdef float var3
    cdef float var4

    for dim_1 in range(dim1):
        for ii in range(1, ns_dim-1):
            for jj in range(1, we_dim-1):
                if np.isnan(grid_mask_view[ii, jj]):
                    var1 = var_out_view[dim_1, ii - 1, jj]
                    var2 = var_out_view[dim_1, ii + 1, jj]
                    var3 = var_out_view[dim_1, ii, jj - 1]
                    var4 = var_out_view[dim_1, ii, jj + 1]
                    var_out_view[dim_1, ii, jj] = 0.25 * (var1 + var2 + var3 + var4)
    # below fixes the "edge" problem caused by destaggering

    return var_out[:, 1:ns_dim-1, 1:we_dim-1]

@cython.boundscheck(False)
@cython.wraparound(False)
def regrid_NMM_rank4(np.ndarray[float, ndim=4] varin, np.ndarray[float, ndim=2] grid_mask, str unstagDim):
    """
    Function to regrid the NMM data when dealing with
    4-dimenional arrays
    """
    begin = datetime.now()
    cdef float[:, :, :, :] varin_view = varin
    cdef int n_dims = len(np.shape(varin))
    cdef int dim_ns_in = np.shape(varin)[n_dims - 2]
    cdef int dim_we_in = np.shape(varin)[n_dims - 1]
    cdef int dim_ns_out = dim_ns_in
    cdef int dim_we_out = dim_we_in * 2
    cdef Py_ssize_t dim1 = np.shape(varin)[0]
    cdef Py_ssize_t dim2 = np.shape(varin)[1]
    cdef Py_ssize_t dim_ns = dim_ns_out
    cdef Py_ssize_t dim_we = dim_we_out

    cdef float[:, :] grid_mask_view = grid_mask

    var_out = np.zeros((dim1, dim2, dim_ns, dim_we), dtype=DTYPE)
    cdef float[:, :, :, :] var_out_view = var_out
    var_out_view[:, :, :, :] = np.nan
    cdef Py_ssize_t ii, jj
    cdef Py_ssize_t ii_count = 0
    cdef Py_ssize_t jj_count = 0
    for ii in range(dim_ns):
        for jj in range(dim_we):
            if not np.isnan(grid_mask_view[ii, jj]):
                var_out_view[:, :, ii, jj] = varin_view[:, :, ii_count, jj_count]
                jj_count = jj_count + 1
                if jj_count == dim_we_in:
                    jj_count = 0
                    ii_count = ii_count + 1

    print(f'total time for fill was {(datetime.now() - begin).total_seconds()}')
    begin = datetime.now()
    # now use your copied values to fill in the blank spots
    #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
    cdef Py_ssize_t ns_dim = np.shape(grid_mask)[0]
    cdef Py_ssize_t we_dim = np.shape(grid_mask)[1]
    cdef Py_ssize_t dim_1

    test = np.zeros((4, dim2), dtype=DTYPE)
    cdef float[:, :] test_view = test
    for dim_1 in range(dim1):
        for ii in range(1, ns_dim-1):
            for jj in range(1, we_dim-1):
                if np.isnan(grid_mask_view[ii, jj]):
                    test_view[0, :] = var_out_view[dim_1, :, ii - 1, jj]
                    test_view[1, :] = var_out_view[dim_1, :, ii + 1, jj]
                    test_view[2, :] = var_out_view[dim_1, :, ii, jj - 1]
                    test_view[3, :] = var_out_view[dim_1, :, ii, jj + 1]
                    var_out[dim_1, :, ii, jj] = 0.25 * np.sum(test_view, 0)
    # below fixes the "edge" problem caused by destaggering

    print(f'total time for interp was {(datetime.now() - begin).total_seconds()}')
    return var_out[:, :, 1:ns_dim-1, 1:we_dim-1]


def wrf_user_unstagger_NMM(varin, unstagDim):
    if unstagDim == "Z":
        var_out = unstagger_NMM_Z(varin)
        # Incase we Destaggered Z## we need to also do the we dim##
        varin = var_out

    n_dims = len(np.shape(varin))
    print(unstagDim, n_dims)

    if n_dims == 2:
        # we only have time and height so no need for spatial fix
        return var_out
    #if n_dims == 2:
    #    # will never hit as is coded now
    #    return regrid_NMM_rank2(varin, unstagDim)

    n_dims = len(np.shape(varin))
    dim_ns_in = np.shape(varin)[n_dims - 2]
    dim_we_in = np.shape(varin)[n_dims - 1]
    dim_ns_out = dim_ns_in
    dim_we_out = dim_we_in * 2
    grid_mask = get_gridmask_NMM(dim_ns_out, dim_we_out, unstagDim)

    if n_dims == 3:
        return regrid_NMM_rank3(varin, grid_mask, unstagDim)
    if n_dims == 4:
        return regrid_NMM_rank4(varin, grid_mask, unstagDim)
