#!/usr/bin/python
import numpy as np


def wrf_user_unstagger_ARW (varin, unstagDim):
    dims =()
    dims = np.shape(varin)
    nd = np.shape(dims)[0]

    if ( unstagDim == "X") | (unstagDim == "U" ):
        dimU = dims[nd-1]
        if ( nd == 5 ):
            varout = 0.5*(varin[:,:,:,:,:dimU-2] + varin[:,:,:,:,1:dimU-0])

        if ( nd == 4 ):
            varout = 0.5*(varin[:,:,:,:dimU-1] + varin[:,:,:,1:dimU-0])

        if ( nd == 3 ):
            varout = 0.5*(varin[:,:,:dimU-1] + varin[:,:,1:dimU-0])

        if ( nd == 2 ):
            varout = 0.5*(varin[:,:dimU-1] + varin[:,1:dimU-0])



    if ( unstagDim == "Y") | (unstagDim == "V" ):
        dimV = dims[nd-2]
        if ( nd == 5 ):
            varout = 0.5*(varin[:,:,:,:dimV-1,:]+varin[:,:,:,1:dimV-0,:])

        if ( nd == 4 ):
            varout = 0.5*(varin[:,:,:dimV-1,:]+varin[:,:,1:dimV-0,:])

        if ( nd == 3 ):
            varout = 0.5*(varin[:,:dimV-1,:]+varin[:,1:dimV-0,:])

        if ( nd == 2 ):
            varout = 0.5*(varin[:dimV-1,:]+varin[1:dimV-0,:])

    if ( unstagDim == "Z" ):
        dimW = dims[nd-3]
        if ( nd == 5 ):
            varout = 0.5*(varin[:,:,0:dimW-1,:,:]+varin[:,:,1:dimW-0,:,:])

        if ( nd == 4 ):
            varout = 0.5*(varin[:,0:dimW-1,:,:]+varin[:,1:dimW-0,:,:])

        if ( nd == 3 ):
            varout = 0.5*(varin[0:dimW-1,:,:]+varin[1:dimW-0,:,:])

    if ( unstagDim == "M" ):
        print('staggared dim is "M", no support for this yet')
        varout = varin


    for ele in ["X","U","Y","V","Z", "M"]:
        if (unstagDim == ele):
            return varout



def wrf_user_unstagger_NMM(varin,unstagDim,test_method=2):
        verbose=False
        dims_in =()
        dims_in = np.shape(varin)
        rank = np.shape(dims_in)[0]
        if verbose:
            print('simon test, ustaggering')
        n_dims=len(np.shape(varin))

        dim_ns_in=np.shape(varin)[n_dims-2]
        dim_we_in=np.shape(varin)[n_dims-1]

        dim_ns_out=dim_ns_in
        dim_we_out=dim_we_in*2

        if verbose:
            print(n_dims)
            print('shape of input',np.shape(varin))
        done_Z=False
        if ( unstagDim == "Z" ):
            if (n_dims == 2):
                var_out=np.zeros((dim_ns_out,dim_we_out),dtype="float")
            if (n_dims == 3):
                dim1=np.shape(varin)[0]
                var_out=np.zeros((dim1,dim_ns_out,dim_we_out),dtype="float")
            if (n_dims == 4):
                dim1=np.shape(varin)[0]
                dim2=np.shape(varin)[1]
                var_out=np.zeros((dim1,dim2,dim_ns_out,dim_we_out),dtype="float")
            dimW = dims_in[rank-3]
            if ( rank == 5 ):
                var_out = 0.5*(varin[:,:,0:dimW-1,:,:]+varin[:,:,1:dimW-0,:,:])
            if ( rank == 4 ):
                var_out = 0.5*(varin[:,0:dimW-1,:,:]+varin[:,1:dimW-0,:,:])
            if ( rank == 3 ):
                var_out = 0.5*(varin[0:dimW-1,:,:]+varin[1:dimW-0,:,:])
            done_Z=True
            #return var_out

        ###Incase we Destaggered Z## we need to also do the we dim##

        if done_Z:
            varin=var_out
            dims_in = np.shape(varin)
            rank = np.shape(dims_in)[0]
            if verbose:
                print('simon test, ustaggering again, WE this time')
            n_dims=len(np.shape(varin))

            dim_ns_in=np.shape(varin)[n_dims-2]
            dim_we_in=np.shape(varin)[n_dims-1]

            dim_ns_out=dim_ns_in
            dim_we_out=dim_we_in*2

        if verbose:
            print(n_dims)
            print('shape of input',np.shape(varin))


        if (n_dims == 2):
            var_out=np.zeros((dim_ns_out,dim_we_out),dtype="float")
        if (n_dims == 3):
            dim1=np.shape(varin)[0]
            var_out=np.zeros((dim1,dim_ns_out,dim_we_out),dtype="float")
        if (n_dims == 4):
            dim1=np.shape(varin)[0]
            dim2=np.shape(varin)[1]
            var_out=np.zeros((dim1,dim2,dim_ns_out,dim_we_out),dtype="float")

        if verbose:
            print('shape of output',np.shape(var_out))

        #general
        grid_mask=np.zeros((dim_ns_out,dim_we_out),dtype="float")
        grid=np.zeros((dim_ns_out,dim_we_out),dtype="float")

        if (unstagDim == 'H') :
            grid_mask[:,:]='NaN'
        if (unstagDim == 'V'):
            grid_mask[:,:]=1
        if (unstagDim == 'Z'):
            grid_mask[:,:]='NaN'


        for ii in range(dim_ns_out):
            for jj in range(dim_we_out):
                if (np.mod(ii,2) == 0) & (np.mod(jj,2) == 0):
                    if (unstagDim=='H') | (unstagDim=='Z'):
                        grid_mask[ii,jj] = 1
                    if (unstagDim=='V'):
                        grid_mask[ii,jj] = 'NaN'

                if (np.mod(ii,2) == 1) & (np.mod(jj,2) == 1):
                    if (unstagDim == 'H') | (unstagDim =='Z'):
                        grid_mask[ii,jj]= 1
                    if (unstagDim=='V'):
                        grid_mask[ii,jj] = 'NaN'

        if (n_dims == 2):
            var_out[:,:]=grid_mask[:,:]
        if (n_dims == 3):
            for dim_1 in range(dim1):
                var_out[dim_1,:,:] = grid_mask[:,:]
        if (n_dims == 4):
            for dim_1 in range(dim1):
                for dim_2 in range(dim2):
                    var_out[dim_1,dim_2,:,:] = grid_mask[:,:]

        if (n_dims == 2):
            ii_count=0
            jj_count=0
            for ii in range(dim_ns_out):
                for jj in range(dim_we_out):
                    if (np.isnan(grid_mask[ii,jj]) == 1):
                        continue
                    else:
                        var_out[ii,jj]=varin[ii_count,jj_count]
                        jj_count=jj_count+1
                        if (jj_count == dim_we_in):
                            jj_count=0
                            ii_count=ii_count+1

        if (n_dims == 3):
            for dim_1 in range(dim1):
                ii_count=0
                jj_count=0
                for ii in range(dim_ns_out):
                    for jj in range(dim_we_out):
                        if (np.isnan(grid_mask[ii,jj]) == 1):
                            continue
                        else:
                            var_out[dim_1,ii,jj]=varin[dim_1,ii_count,jj_count]
                            jj_count=jj_count+1
                            if (jj_count == dim_we_in):
                                jj_count=0
                                ii_count=ii_count+1


        if (n_dims == 4):
            for dim_1 in range(dim1):
                #for dim_2 in range(dim2):
                if 1==1:
                    ii_count=0
                    jj_count=0
                    for ii in range(dim_ns_out):
                        for jj in range(dim_we_out):
                            if (np.isnan(grid_mask[ii,jj]) == 1):
                                continue
                            else:
                                var_out[dim_1,:,ii,jj]=varin[dim_1,:,ii_count,jj_count]
                                jj_count=jj_count+1
                                if (jj_count == dim_we_in):
                                    jj_count=0
                                    ii_count=ii_count+1

        ###################now use your copied values to fill in the blank spots###########

        if (test_method==1):
            if (n_dims == 2):
                for ii in range(np.shape(grid_mask)[0]):
                    for jj in range(np.shape(grid_mask)[1]):
                        if (np.isnan(var_out[ii,jj]) == 1):
                            if (jj+1 == np.shape(var_out)[n_dims-1]) | (jj-1 ==-1): 
                                continue
                            else:
                                var_out[ii,jj]=0.5*(var_out[ii,jj-1]+var_out[ii,jj+1])

            if (n_dims == 3):
                for dim_1 in range(dim1):
                    for ii in range(np.shape(grid_mask)[0]):
                        for jj in range(np.shape(grid_mask)[1]):
                            if (np.isnan(var_out[dim_1,ii,jj]) == 1):
                                if (jj+1 == np.shape(var_out)[n_dims-1]) | (jj-1 ==-1): 
                                    continue
                                else:
                                    var_out[dim_1,ii,jj]=0.5*(var_out[dim_1,ii,jj-1]+var_out[dim_1,ii,jj+1])

            if (n_dims == 4):
                for dim_1 in range(dim1):
                    for dim_2 in range(dim2):
                        for ii in range(np.shape(grid_mask)[0]):
                            for jj in range(np.shape(grid_mask)[1]):
                                if (np.isnan(var_out[dim_1,dim_2,ii,jj]) == 1):
                                    if (jj+1 == np.shape(var_out)[n_dims-1]) | (jj-1 ==-1): 
                                        continue
                                    else:
                                        var_out[dim_1,dim_2,ii,jj]=0.5*(var_out[dim_1,dim_2,ii,jj-1]+var_out[dim_1,dim_2,ii,jj+1])

        
        if (test_method==2):
#        try:
            #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
            if (n_dims == 2):
                for ii in range(np.shape(grid_mask)[0]):
                    for jj in range(np.shape(grid_mask)[1]):
                        if (np.isnan(var_out[ii,jj])):
                            if (ii+1 == np.shape(var_out)[0]) | (ii-1 ==-1): 
                                continue
                            elif (jj+1 == np.shape(var_out)[1]) | (jj-1 ==-1): 
                                continue
                            else:
                                var_out[ii,jj]=0.25*(var_out[ii-1,jj]+var_out[ii+1,jj]+var_out[ii,jj-1]+var_out[ii,jj+1])
            if (n_dims == 3):
                for ii in range(np.shape(grid_mask)[0]):
                    for jj in range(np.shape(grid_mask)[1]):
                        if (np.isnan(var_out[0,ii,jj])):
                            if 1==1:
                                try:
                                    var1=var_out[0,ii-1,jj]
                                except:
                                    var1=np.nan
                                try:
                                    var2=var_out[0,ii+1,jj]
                                except:
                                    var2=np.nan
                                try:
                                    var3=var_out[0,ii,jj-1]
                                except:
                                    var3=np.nan
                                try:
                                    var4=var_out[0,ii,jj+1]
                                except:
                                    var4=np.nan
                                var_list=[var1,var2,var3,var4]
                                nan_sum=np.nansum(var_list)
                                valid_count=(len(var_list)-np.sum(np.isnan(var_list)))
                                var_out[0,ii,jj]=nan_sum/valid_count

            if (n_dims == 4):
                test=np.zeros((4,dim2))
                for dim_1 in range(dim1):
                    #for dim_2 in range(dim2):
                        for ii in range(np.shape(grid_mask)[0]):
                            if verbose:
                                print(ii)
                            for jj in range(np.shape(grid_mask)[1]):
                                if (np.isnan(grid_mask[ii,jj])):
                                        try:
                                            test[0,:]=var_out[0,:,ii-1,jj]
                                        except:
                                            test[0,:]=[np.nan for x in range(dim2)]
                                        try:
                                            test[1,:]=var_out[0,:,ii+1,jj]
                                        except:
                                            test[1,:]=[np.nan for x in range(dim2)]
                                        try:
                                            test[2,:]=var_out[0,:,ii,jj-1]
                                        except:
                                            test[2,:]=[np.nan for x in range(dim2)]
                                        try:
                                            test[3,:]=var_out[0,:,ii,jj+1]
                                        except:
                                            test[3,:]=[np.nan for x in range(dim2)]
                                        nan_sum=np.nansum(test,0)
                                        nan_count=np.sum(np.isnan(test),0)
                                        valid_count=np.shape(test)[0]-nan_count
                                        var_out[0,:,ii,jj]=nan_sum/valid_count
        # below fixes the "edge" problem caused by destaggering
        if (n_dims == 2):
            return var_out[1:-1, 1:-1]
        if (n_dims == 3):
            return var_out[:, 1:-1, 1:-1]
        if (n_dims == 4):
            return var_out[:, :, 1:-1, 1:-1]
 

