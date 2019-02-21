#!/usr/bin/python
import numpy as np


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


def get_gridmask_NMM(dim_ns_out, dim_we_out, unstagDim):
    """
    This part is general and should not depend on the
    higher dimenions, essentially setting up how the 
    grid is structured
    """

    grid_mask = np.zeros((dim_ns_out, dim_we_out), dtype="float")

    if (unstagDim == 'H'):
        grid_mask[:, :] = np.nan
    if (unstagDim == 'V'):
        grid_mask[:, :] = 1
    if (unstagDim == 'Z'):
        grid_mask[:, :] = np.nan

    for ii in range(dim_ns_out):
        for jj in range(dim_we_out):
            if np.mod(ii, 2) == 0 and np.mod(jj, 2) == 0:
                if unstagDim == 'H' or unstagDim == 'Z':
                    grid_mask[ii, jj] = 1
                elif unstagDim == 'V':
                    grid_mask[ii, jj] = np.nan

            if np.mod(ii, 2) == 1 and np.mod(jj, 2) == 1:
                if unstagDim == 'H' or unstagDim == 'Z':
                    grid_mask[ii, jj] = 1
                if unstagDim == 'V':
                    grid_mask[ii, jj] = np.nan
    return grid_mask


def regrid_NMM_rank2(varin, unstagDim):
    """
    Function to regrid the NMM data when dealing with
    2-dimenional arrays.
    Actually this probably never occurs due to the fact that if
    we have 2-d data one of them is probably time or level
    """

    n_dims = len(np.shape(varin))
    # spatial dimensions are always the last two
    dim_ns_in = np.shape(varin)[n_dims - 2]
    dim_we_in = np.shape(varin)[n_dims - 1]
    # output grid will have double the WE dimension
    dim_ns_out = dim_ns_in
    dim_we_out = dim_we_in * 2

    grid_mask = get_gridmask_NMM(dim_ns_out, dim_we_out, unstagDim)

    # 2 times but no time, i dont think this ever occurs
    var_out = np.zeros((dim_ns_out, dim_we_out), dtype="float")
    var_out[:, :] = grid_mask[:, :]
    ii_count = 0
    jj_count = 0
    for ii in range(dim_ns_out):
        for jj in range(dim_we_out):
            if not np.isnan(grid_mask[ii, jj]):
                var_out[ii, jj] = varin[ii_count, jj_count]
                jj_count = jj_count + 1
                if jj_count == dim_we_in:
                    jj_count = 0
                    ii_count = ii_count + 1
    # now use your copied values to fill in the blank spots
    #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
    ns_dim = np.shape(grid_mask)[0]
    we_dim = np.shape(grid_mask)[1]
    for ii in range(1,  ns_dim-1):
        for jj in range(1, we_dim-1):
            if np.isnan(grid_mask[ii, jj]):
                var1 = var_out[ii - 1, jj]
                var2 = var_out[ii + 1, jj]
                var3 = var_out[ii, jj - 1]
                var4 = var_out[ii, jj + 1]
                var_out[ii, jj] = 0.25 * (var1 + var2 + var3 + var4)
    # below fixes the "edge" problem caused by destaggering
    return var_out[1:ns_dim-1, 1:we_dim-1]

def regrid_NMM_rank3(varin, unstagDim):
    """
    Function to regrid the NMM data when dealing with
    3-dimenional arrays
    """

    n_dims = len(np.shape(varin))
    # spatial dimensions are always the last two
    dim_ns_in = np.shape(varin)[n_dims - 2]
    dim_we_in = np.shape(varin)[n_dims - 1]
    # output grid will have double the WE dimension
    dim_ns_out = dim_ns_in
    dim_we_out = dim_we_in * 2

    grid_mask = get_gridmask_NMM(dim_ns_out, dim_we_out, unstagDim)

    ##
    dim1 = np.shape(varin)[0]
    var_out = np.zeros((dim1, dim_ns_out, dim_we_out), dtype="float")
    var_out[:, :, :] = np.nan
    ii_count = 0
    jj_count = 0
    for ii in range(dim_ns_out):
        for jj in range(dim_we_out):
            if not np.isnan(grid_mask[ii, jj]):
                var_out[:, ii, jj] = varin[:, ii_count, jj_count]
                jj_count = jj_count + 1
                if (jj_count == dim_we_in):
                    jj_count = 0
                    ii_count = ii_count + 1
    # now use your copied values to fill in the blank spots
    #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
    ns_dim = np.shape(grid_mask)[0]
    we_dim = np.shape(grid_mask)[1]
    for dim_1 in range(dim1):
        for ii in range(1, ns_dim-1):
            for jj in range(1, we_dim-1):
                if np.isnan(grid_mask[ii, jj]):
                    var1 = var_out[dim_1, ii - 1, jj]
                    var2 = var_out[dim_1, ii + 1, jj]
                    var3 = var_out[dim_1, ii, jj - 1]
                    var4 = var_out[dim_1, ii, jj + 1]
                    var_out[dim_1, ii, jj] = 0.25 * (var1 + var2 + var3 + var4)
    # below fixes the "edge" problem caused by destaggering
    return var_out[:, 1:ns_dim-1, 1:we_dim-1]

def regrid_NMM_rank4(varin, unstagDim):
    """
    Function to regrid the NMM data when dealing with
    4-dimenional arrays
    """

    n_dims = len(np.shape(varin))
    # spatial dimensions are always the last two
    dim_ns_in = np.shape(varin)[n_dims - 2]
    dim_we_in = np.shape(varin)[n_dims - 1]
    # output grid will have double the WE dimension
    dim_ns_out = dim_ns_in
    dim_we_out = dim_we_in * 2

    grid_mask = get_gridmask_NMM(dim_ns_out, dim_we_out, unstagDim)

    ##
    dim1 = np.shape(varin)[0]
    dim2 = np.shape(varin)[1]
    var_out = np.zeros((dim1, dim2, dim_ns_out, dim_we_out), dtype="float")
    var_out[:, :, :, :] = np.nan
    ii_count = 0
    jj_count = 0
    for ii in range(dim_ns_out):
        for jj in range(dim_we_out):
            if not np.isnan(grid_mask[ii, jj]):
                var_out[:, :, ii, jj] = varin[:, :, ii_count, jj_count]
                jj_count = jj_count + 1
                if jj_count == dim_we_in:
                    jj_count = 0
                    ii_count = ii_count + 1
    # now use your copied values to fill in the blank spots
    #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
    ns_dim = np.shape(grid_mask)[0]
    we_dim = np.shape(grid_mask)[1]
    test = np.zeros((4, dim2))
    for dim_1 in range(dim1):
        for ii in range(1, ns_dim-1):
            for jj in range(1, we_dim-1):
                if np.isnan(grid_mask[ii, jj]):
                    test[0, :] = var_out[dim_1, :, ii - 1, jj]
                    test[1, :] = var_out[dim_1, :, ii + 1, jj]
                    test[2, :] = var_out[dim_1, :, ii, jj - 1]
                    test[3, :] = var_out[dim_1, :, ii, jj + 1]
                    var_out[dim_1, :, ii, jj] = 0.25 * np.sum(test, 0)
    # below fixes the "edge" problem caused by destaggering
    return var_out[:, :, 1:ns_dim-1, 1:we_dim-1]


def wrf_user_unstagger_NMM(varin, unstagDim):
    done_Z = False
    if unstagDim == "Z":
        var_out = unstagger_NMM_Z(varin)
        done_Z = True
        # Incase we Destaggered Z## we need to also do the we dim##
        varin = var_out

    n_dims = len(np.shape(varin))
    print(unstagDim, n_dims)
    if n_dims == 2:
        # we only have time and height so no need for spatial fix
        return var_out
    if n_dims == 2:
        # will never hit as is coded now
        return regrid_NMM_rank2(varin, unstagDim)
    if n_dims == 3:
        return regrid_NMM_rank3(varin, unstagDim)
    if n_dims == 4:
        return regrid_NMM_rank4(varin, unstagDim)
