#!/usr/bin/python
import numpy as np
cimport numpy as np
cimport cython
#from datetime import datetime

DTYPE = np.float32

def wrf_user_unstagger_ARW(varin, unstagDim):
    #begin = datetime.now()
    dims = np.shape(varin)
    nd = np.shape(dims)[0]

    #print(unstagDim, nd)
    if unstagDim == "X" or unstagDim == "U":
        dimU = dims[nd - 1]
        if nd == 5:
            varout = 0.5 * (
                varin[:, :, :, :, :dimU - 2] + varin[:, :, :, :, 1:dimU - 0])
        elif nd == 4:
            varout = 0.5 * (
                varin[:, :, :, :dimU - 1] + varin[:, :, :, 1:dimU - 0])
        elif nd == 3:
            varout = 0.5 * (varin[:, :, :dimU - 1] + varin[:, :, 1:dimU - 0])
        elif nd == 2:
            varout = 0.5 * (varin[:, :dimU - 1] + varin[:, 1:dimU - 0])
    elif unstagDim == "Y" or unstagDim == "V":
        dimV = dims[nd - 2]
        if nd == 5:
            varout = 0.5 * (
                varin[:, :, :, :dimV - 1, :] + varin[:, :, :, 1:dimV - 0, :])
        elif nd == 4:
            varout = 0.5 * (
                varin[:, :, :dimV - 1, :] + varin[:, :, 1:dimV - 0, :])
        elif nd == 3:
            varout = 0.5 * (varin[:, :dimV - 1, :] + varin[:, 1:dimV - 0, :])
        elif nd == 2:
            varout = 0.5 * (varin[:dimV - 1, :] + varin[1:dimV - 0, :])
    elif unstagDim == "Z":
        dimW = dims[nd - 3]
        if nd == 5:
            varout = 0.5 * (
                varin[:, :, 0:dimW - 1, :, :] + varin[:, :, 1:dimW - 0, :, :])
        elif nd == 4:
            varout = 0.5 * (
                varin[:, 0:dimW - 1, :, :] + varin[:, 1:dimW - 0, :, :])
        elif nd == 3:
            varout = 0.5 * (varin[0:dimW - 1, :, :] + varin[1:dimW - 0, :, :])
        elif nd == 2:
            varout = 0.5 * (varin[:, 0:dimW - 1] + varin[:, 1:dimW - 0])
    elif unstagDim == "M":
        print('staggared dim is "M", no support for this yet')
        varout = varin

    #print(f'total time for ARW destagger was {(datetime.now() - begin).total_seconds()}')
    for ele in ["X", "U", "Y", "V", "Z", "M"]:
        if (unstagDim == ele):
            return varout

def unstagger_NMM_Z(varin):
    #begin = datetime.now()
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
    #print(f'total time for unstag_Z was {(datetime.now() - begin).total_seconds()}')
    return var_out


@cython.boundscheck(False)
@cython.wraparound(False)
def get_gridmask_NMM(Py_ssize_t dim_ns, Py_ssize_t dim_we, str unstagDim):
    """
    This part is general and should not depend on the
    higher dimenions, essentially setting up how the 
    grid is structured
    """

    #begin = datetime.now()
    grid_mask = np.zeros((dim_ns, dim_we), dtype=DTYPE)
    cdef float[:, :] grid_mask_view = grid_mask
    cdef float NaN = float('NaN')
    cdef Py_ssize_t ii, jj

    if unstagDim == 'H' or unstagDim == 'Z':
        grid_mask_view[:, :] = NaN
    elif unstagDim == 'V':
        grid_mask_view[:, :] = 1

    if unstagDim == 'H' or unstagDim == 'Z':
        for ii in range(dim_ns):
            for jj in range(dim_we):
                if ii % 2 == 0 and jj %  2 == 0:
                    grid_mask_view[ii, jj] = 1
                elif ii % 2 == 1 and jj % 2 == 1:
                    grid_mask_view[ii, jj] = 1
    elif unstagDim == 'V':
        for ii in range(dim_ns):
            for jj in range(dim_we):
                if ii % 2 == 0 and jj % 2 == 0:
                    grid_mask_view[ii, jj] = NaN
                elif ii % 2 == 1 and jj % 2 == 1:
                    grid_mask_view[ii, jj] = NaN
    #print(f'total time for grid_mask was {(datetime.now() - begin).total_seconds()}')
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

@cython.boundscheck(False)
@cython.wraparound(False)
def regrid_NMM_rank3(np.ndarray[float, ndim=3] varin, np.ndarray[float, ndim=2] grid_mask, str unstagDim):
    """
    Function to regrid the NMM data when dealing with
    3-dimenional arrays
    """

    #begin = datetime.now()
    cdef float NaN = float('NaN')
    cdef int dim_we_orig = np.shape(varin)[2]
    cdef Py_ssize_t dim1 = np.shape(varin)[0] # typically time but not neccessarly
    cdef Py_ssize_t dim_ns = np.shape(varin)[1]
    cdef Py_ssize_t dim_we = np.shape(varin)[2] * 2
    cdef Py_ssize_t indx_1 # dim1 index
    cdef Py_ssize_t ii, jj # spatial index
    cdef Py_ssize_t ii_count = 0 # spatial indx on orig grid
    cdef Py_ssize_t jj_count = 0 # spatial indx on orig grid
    # defined view_arrays for fast access
    cdef float[:, :, :] varin_view = varin
    cdef float[:, :] grid_mask_view = grid_mask
    var_out = np.zeros((dim1, dim_ns, dim_we), dtype=DTYPE)
    cdef float[:, :, :] var_out_view = var_out
    var_out_view[:, :, :] = NaN

    for indx_1 in range(dim1):
        ii_count = 0
        jj_count = 0
        for ii in range(dim_ns):
            for jj in range(dim_we):
                # nans and equality are tricky if we are not using is.nan
                if grid_mask_view[ii, jj] == 1:
                    var_out_view[indx_1, ii, jj] = varin_view[indx_1, ii_count, jj_count]
                    jj_count = jj_count + 1
                    if jj_count == dim_we_orig:
                        jj_count = 0
                        ii_count = ii_count + 1

    #print(f'total time for fill was {(datetime.now() - begin).total_seconds()}')
    #begin = datetime.now()
    # now use your copied values to fill in the blank spots
    #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
    cdef float var1
    cdef float var2
    cdef float var3
    cdef float var4

    for indx_1 in range(dim1):
        for ii in range(1, dim_ns-1):
            for jj in range(1, dim_we-1):
                if grid_mask_view[ii, jj] != 1:
                    var1 = var_out_view[indx_1, ii - 1, jj]
                    var2 = var_out_view[indx_1, ii + 1, jj]
                    var3 = var_out_view[indx_1, ii, jj - 1]
                    var4 = var_out_view[indx_1, ii, jj + 1]
                    var_out_view[indx_1, ii, jj] = 0.25 * (var1 + var2 + var3 + var4)
    # below fixes the "edge" problem caused by destaggering
    #print(f'total time for interp was {(datetime.now() - begin).total_seconds()}')
    return var_out[:, 1:dim_ns-1, 1:dim_we-1]

@cython.boundscheck(False)
@cython.wraparound(False)
def regrid_NMM_rank4(np.ndarray[float, ndim=4] varin, np.ndarray[float, ndim=2] grid_mask, str unstagDim):
    """
    Function to regrid the NMM data when dealing with
    4-dimenional arrays
    """
    #begin = datetime.now()
    cdef float NaN = float('NaN')
    cdef Py_ssize_t dim1 = np.shape(varin)[0] # typically time but not neccessarly
    cdef Py_ssize_t dim2 = np.shape(varin)[1] # typically bottom-top but not necessarly
    cdef Py_ssize_t dim_ns = np.shape(varin)[2]
    cdef int dim_we_orig = np.shape(varin)[3]
    cdef Py_ssize_t dim_we = np.shape(varin)[3] * 2
    cdef Py_ssize_t indx_1 # dim1 index
    cdef Py_ssize_t indx_2 # dim2 index
    cdef Py_ssize_t ii, jj # spatial index
    cdef Py_ssize_t ii_count = 0 # spatial indx on orig grid
    cdef Py_ssize_t jj_count = 0 # spatial indx on orig grid

    # defined view_arrays for fast access
    cdef float[:, :, :, :] varin_view = varin
    cdef float[:, :] grid_mask_view = grid_mask
    var_out = np.zeros((dim1, dim2, dim_ns, dim_we), dtype=DTYPE)
    cdef float[:, :, :, :] var_out_view = var_out
    var_out_view[:, :, :, :] = NaN

    for indx_1 in range(dim1):
        for indx_2 in range(dim2):
            ii_count = 0
            jj_count = 0
            for ii in range(dim_ns):
                for jj in range(dim_we):
                    if grid_mask_view[ii, jj] == 1:
                        var_out_view[indx_1, indx_2, ii, jj] = varin_view[indx_1, indx_2, ii_count, jj_count]
                        jj_count = jj_count + 1
                        if jj_count == dim_we_orig:
                            jj_count = 0
                            ii_count = ii_count + 1

    #print(f'total time for fill was {(datetime.now() - begin).total_seconds()}')
    #begin = datetime.now()
    # now use your copied values to fill in the blank spots
    #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
    cdef float var1
    cdef float var2
    cdef float var3
    cdef float var4

    for indx_1 in range(dim1):
        for indx_2 in range(dim2):
            for ii in range(1, dim_ns-1):
                for jj in range(1, dim_we-1):
                    if grid_mask_view[ii, jj] != 1:
                        var1 = var_out_view[indx_1, indx_2, ii - 1, jj]
                        var2 = var_out_view[indx_1, indx_2, ii + 1, jj]
                        var3 = var_out_view[indx_1, indx_2, ii, jj - 1]
                        var4 = var_out_view[indx_1, indx_2, ii, jj + 1]
                        var_out_view[indx_1, indx_2, ii, jj] = 0.25 * (var1 + var2 + var3 + var4)
    # below fixes the "edge" problem caused by destaggering

    #print(f'total time for interp was {(datetime.now() - begin).total_seconds()}')
    return var_out[:, :, 1:dim_ns-1, 1:dim_we-1]


def wrf_user_unstagger_NMM(varin, unstagDim):
    if unstagDim == "Z":
        var_out = unstagger_NMM_Z(varin)
        # Incase we Destaggered Z## we need to also do the we dim##
        varin = var_out

    n_dims = len(np.shape(varin))
    #print(unstagDim, n_dims)

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
