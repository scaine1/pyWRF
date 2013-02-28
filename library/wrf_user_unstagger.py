#!/usr/bin/python
import numpy as n


def wrf_user_unstagger_ARW (varin, unstagDim):
    dims =()
    dims = n.shape(varin)
    nd = n.shape(dims)[0]

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
	print 'staggared dim is "M", no support for this yet'
	varout = varin


    for ele in ["X","U","Y","V","Z", "M"]:
	if (unstagDim == ele):
	    return varout



def wrf_user_unstagger_NMM(varin,unstagDim,test_method=1):
	dims =()
	dims = n.shape(varin)
	nd = n.shape(dims)[0]
	print 'simon test'
	print nd

	n_dims=len(n.shape(varin))
	dim_ii=n.shape(varin)[n_dims-2]
	dim_jj=n.shape(varin)[n_dims-1]
	
	ii_dim_h=dim_ii
	jj_dim_h=dim_jj*2

	if ( unstagDim == "Z" ):

	    if (n_dims == 2):
		var_out=n.zeros((dim_ii,dim_jj),dtype="float")

	    if (n_dims == 3):
		dim1=n.shape(varin)[0]
		var_out=n.zeros((dim1,dim_ii,dim_jj),dtype="float")

	    if (n_dims == 4):
		dim1=n.shape(varin)[0]
		dim2=n.shape(varin)[1]
		var_out=n.zeros((dim1,dim2,dim_ii,dim_jj),dtype="float")


	    dimW = dims[nd-3]
	    if ( nd == 5 ):
		var_out = 0.5*(varin[:,:,0:dimW-1,:,:]+varin[:,:,1:dimW-0,:,:])

	    if ( nd == 4 ):
		var_out = 0.5*(varin[:,0:dimW-1,:,:]+varin[:,1:dimW-0,:,:])

	    if ( nd == 3 ):
		var_out = 0.5*(varin[0:dimW-1,:,:]+varin[1:dimW-0,:,:])

	    return var_out



	
	if (n_dims == 2):
	    var_out=n.zeros((ii_dim_h,jj_dim_h),dtype="float")

	if (n_dims == 3):
	    dim1=n.shape(varin)[0]
	    var_out=n.zeros((dim1,ii_dim_h,jj_dim_h),dtype="float")

	if (n_dims == 4):
	    dim1=n.shape(varin)[0]
	    dim2=n.shape(varin)[1]
	    var_out=n.zeros((dim1,dim2,ii_dim_h,jj_dim_h),dtype="float")





	#general
	grid_mask=n.zeros((ii_dim_h,jj_dim_h),dtype="float")
	grid=n.zeros((ii_dim_h,jj_dim_h),dtype="float")
	
	
	if (unstagDim == 'H'):
	    grid_mask[:,:]='NaN'
	if (unstagDim == 'V'):
	    grid_mask[:,:]=1
	
	for ii in range(ii_dim_h):
	    for jj in range(jj_dim_h):
		if (n.mod(ii,2) == 0) & (n.mod(jj,2) == 0):
		    if (unstagDim=='H'):
			grid_mask[ii,jj] = 1
		    if (unstagDim=='V'):
			grid_mask[ii,jj] = 'NaN'
	
		if (n.mod(ii,2) == 1) & (n.mod(jj,2) == 1):
		    if (unstagDim == 'H'):
			grid_mask[ii,jj]= 1
		    if (unstagDim=='V'):
			grid_mask[ii,jj] = 'NaN'
	
	

	print n_dims
	print n.shape(varin)


	if (n_dims == 2):
	    var_out[:,:]=grid_mask[:,:]
	
	if (n_dims == 3):
	    for dim_1 in range(dim1):
		var_out[dim_1,:,:] = grid_mask[:,:]

	if (n_dims == 4):
	    for dim_1 in range(dim1):
		for dim_2 in range(dim2):
		    var_out[dim_1,dim_2,:,:] = grid_mask[:,:]

	
	ii_count=0
	jj_count=0



	if (n_dims == 2):
	    for ii in range(ii_dim_h):
		for jj in range(jj_dim_h):
		    if (n.isnan(grid_mask[ii,jj]) == 1):
			continue
		    else:
			var_out[ii,jj]=varin[ii_count,jj_count]
			jj_count=jj_count+1
			if (jj_count == dim_jj):
			    jj_count=0
			    ii_count=ii_count+1
	

	if (n_dims == 3):
	    for dim_1 in range(dim1):
		for ii in range(ii_dim_h):
		    for jj in range(jj_dim_h):
			if (n.isnan(grid_mask[ii,jj]) == 1):
			    continue
			else:
			    var_out[dim_1,ii,jj]=varin[dim_1,ii_count,jj_count]
			    jj_count=jj_count+1
			    if (jj_count == dim_jj):
				jj_count=0
				ii_count=ii_count+1


	if (n_dims == 4):
	    for dim_1 in range(dim1):
		for dim_2 in range(dim2):
		    for ii in range(ii_dim_h):
			for jj in range(jj_dim_h):
			    if (n.isnan(grid_mask[ii,jj]) == 1):
				continue
			    else:
				var_out[dim_1,dim_2,ii,jj]=varin[dim_1,dim_2,ii_count,jj_count]
				jj_count=jj_count+1
				if (jj_count == dim_jj):
				    jj_count=0
				    ii_count=ii_count+1



	
	###################now use your copied values to fill in the blank spots###########
	
	
	if (test_method==1):
	    if (n_dims == 2):
		for ii in range(n.shape(grid_mask)[0]):
		    for jj in range(n.shape(grid_mask)[1]):
			if (n.isnan(var_out[ii,jj]) == 1):
			    if (jj+1 == n.shape(var_out)[n_dims-1]) | (jj-1 ==-1): 
				continue
			    else:
				var_out[ii,jj]=0.5*(var_out[ii,jj-1]+var_out[ii,jj+1])
	
	    if (n_dims == 3):
		for dim_1 in range(dim1):
		    for ii in range(n.shape(grid_mask)[0]):
			for jj in range(n.shape(grid_mask)[1]):
			    if (n.isnan(var_out[dim_1,ii,jj]) == 1):
				if (jj+1 == n.shape(var_out)[n_dims-1]) | (jj-1 ==-1): 
				    continue
				else:
				    var_out[dim_1,ii,jj]=0.5*(var_out[dim_1,ii,jj-1]+var_out[dim_1,ii,jj+1])

	    if (n_dims == 4):
		for dim_1 in range(dim1):
		    for dim_2 in range(dim2):
			for ii in range(n.shape(grid_mask)[0]):
			    for jj in range(n.shape(grid_mask)[1]):
				if (n.isnan(var_out[dim_1,dim_2,ii,jj]) == 1):
				    if (jj+1 == n.shape(var_out)[n_dims-1]) | (jj-1 ==-1): 
					continue
				    else:
					var_out[dim_1,dim_2,ii,jj]=0.5*(var_out[dim_1,dim_2,ii,jj-1]+var_out[dim_1,dim_2,ii,jj+1])

	
	if (test_method==2):
	    #this method uses 4 grid ;points to work out the value of the point, the borders have nans though
	    for ii in range(n.shape(grid_mask)[0]):
		for jj in range(n.shape(grid_mask)[1]):
		    if (n.isnan(var_out[ii,jj]) == 1):
			if (ii+1 == n.shape(var_out)[0]) | (ii-1 ==-1): 
			    continue
			elif (jj+1 == n.shape(var_out)[1]) | (jj-1 ==-1): 
			    continue
		        else:
			    var_out[ii,jj]=0.25*(var_out[ii-1,jj]+var_out[ii+1,jj]+var_out[ii,jj-1]+var_out[ii,jj+1])


	return var_out
 

