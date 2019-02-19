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


    #general
def get_gridmask_NMM(dim_ns_out, dim_we_out, unstagDim):
    """
    This part is general and should not depend on the
    higher dimenions, essentially setting up how the 
    grid is structured
    """

    grid_mask = np.zeros((dim_ns_out, dim_we_out), dtype="float")
    grid = np.zeros((dim_ns_out, dim_we_out), dtype="float")

    if (unstagDim == 'H'):
        grid_mask[:, :] = 'NaN'
    if (unstagDim == 'V'):
        grid_mask[:, :] = 1
    if (unstagDim == 'Z'):
        grid_mask[:, :] = 'NaN'

    for ii in range(dim_ns_out):
        for jj in range(dim_we_out):
            if (np.mod(ii, 2) == 0) & (np.mod(jj, 2) == 0):
                if (unstagDim == 'H') | (unstagDim == 'Z'):
                    grid_mask[ii, jj] = 1
                if (unstagDim == 'V'):
                    grid_mask[ii, jj] = 'NaN'

            if (np.mod(ii, 2) == 1) & (np.mod(jj, 2) == 1):
                if (unstagDim == 'H') | (unstagDim == 'Z'):
                    grid_mask[ii, jj] = 1
                if (unstagDim == 'V'):
                    grid_mask[ii, jj] = 'NaN'
    return grid_mask

def wrf_user_unstagger_NMM(varin, unstagDim):
    done_Z = False
    if (unstagDim == "Z"):
        var_out = unstagger_NMM_Z(varin)
        done_Z = True

    ###Incase we Destaggered Z## we need to also do the we dim##
    if done_Z:
        varin = var_out

    dims_in = np.shape(varin)
    n_dims = len(dims_in)
    print(unstagDim, n_dims)
    if n_dims == 2:
        # we only have time and height so no need for spatial
        # fix
        return var_out

    # spatial dimensions are always the last two
    dim_ns_in = np.shape(varin)[n_dims - 2]
    dim_we_in = np.shape(varin)[n_dims - 1]
    # output grid will have double the WE dimension
    dim_ns_out = dim_ns_in
    dim_we_out = dim_we_in * 2

    grid_mask = get_gridmask_NMM(dim_ns_out, dim_we_out, unstagDim)

    if (n_dims == 2):
        # 2 times but no time, i dont think this ever occurs
        var_out = np.zeros((dim_ns_out, dim_we_out), dtype="float")
        var_out[:, :] = grid_mask[:, :]
        ii_count = 0
        jj_count = 0
        for ii in range(dim_ns_out):
            for jj in range(dim_we_out):
                if np.isnan(grid_mask[ii, jj]):
                    continue
                else:
                    var_out[ii, jj] = varin[ii_count, jj_count]
                    jj_count = jj_count + 1
                    if jj_count == dim_we_in:
                        jj_count = 0
                        ii_count = ii_count + 1
        # now use your copied values to fill in the blank spots
        #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
        for ii in range(np.shape(grid_mask)[0]):
            for jj in range(np.shape(grid_mask)[1]):
                #if np.isnan(var_out[ii, jj]):
                if np.isnan(grid_mask[ii, jj]):
                    if (ii + 1 == np.shape(var_out)[0]) | (ii - 1 == -1):
                        continue
                    elif (jj + 1 == np.shape(var_out)[1]) | (jj - 1 == -1):
                        continue
                    else:
                        var_out[ii, jj] = 0.25 * (
                            var_out[ii - 1, jj] + var_out[ii + 1, jj] +
                            var_out[ii, jj - 1] + var_out[ii, jj + 1])


    if (n_dims == 3):
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
        for dim_1 in range(dim1):
            for ii in range(1, np.shape(grid_mask)[0] -1):
                for jj in range(1, np.shape(grid_mask)[1] -1):
                    if np.isnan(grid_mask[ii, jj]):
                        var1 = var_out[dim_1, ii - 1, jj]
                        var2 = var_out[dim_1, ii + 1, jj]
                        var3 = var_out[dim_1, ii, jj - 1]
                        var4 = var_out[dim_1, ii, jj + 1]
                        var_out[dim_1, ii, jj] = (var1 + var2 + var3 + var4) /  4.

    if (n_dims == 4):
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
        test = np.zeros((4, dim2))
        for dim_1 in range(dim1):
            for ii in range(1, np.shape(grid_mask)[0] -1):
                for jj in range(1, np.shape(grid_mask)[1] -1):
                    if np.isnan(grid_mask[ii, jj]):
                        test[0, :] = var_out[dim_1, :, ii - 1, jj]
                        test[1, :] = var_out[dim_1, :, ii + 1, jj]
                        test[2, :] = var_out[dim_1, :, ii, jj - 1]
                        test[3, :] = var_out[dim_1, :, ii, jj + 1]
                        var_out[dim_1, :, ii, jj] = np.sum(test, 0) / 4.

    # Having memory issues with big runs
    # testing reseting of variables for Garbage Collection
    # not sure if this actually does anything
    test = None
    grid_mask = None
    grid = None
    varin = None

    # below fixes the "edge" problem caused by destaggering
    if (n_dims == 2):
        return var_out[1:-1, 1:-1]
    if (n_dims == 3):
        return var_out[:, 1:-1, 1:-1]
    if (n_dims == 4):
        return var_out[:, :, 1:-1, 1:-1]
