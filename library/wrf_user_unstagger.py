#!/usr/bin/python
import numpy as n


def wrf_user_unstagger (varin, unstagDim):
    dims =()
    dims = n.shape(varin)
    nd = n.shape(dims)[0]

    if ( unstagDim == "X") | (unstagDim == "U" ):
	dimU = dims[nd-1]
	if ( nd == 5 ):
	    #varout = 0.5*(varin[:,:,:,:,:dimU-2] + varin[:,:,:,:,1:dimU-1])  #python indexing appears to work differently than ncl, 1 less
	    varout = 0.5*(varin[:,:,:,:,:dimU-2] + varin[:,:,:,:,1:dimU-0])

	if ( nd == 4 ):
	    #varout = 0.5*(varin[:,:,:,:dimU-2] + varin[:,:,:,1:dimU-1])
	    varout = 0.5*(varin[:,:,:,:dimU-1] + varin[:,:,:,1:dimU-0])

	if ( nd == 3 ):
	    #varout = 0.5*(varin[:,:,:dimU-2] + varin[:,:,1:dimU-1])
	    varout = 0.5*(varin[:,:,:dimU-1] + varin[:,:,1:dimU-0])

	if ( nd == 2 ):
	    #varout = 0.5*(varin[:,:dimU-2] + varin[:,1:dimU-1])
	    varout = 0.5*(varin[:,:dimU-1] + varin[:,1:dimU-0])

#	do i = 0,nd-2
#	    varout!i = varin!i
#	end do
#	i = nd-1
#	varout!i = "west_east"
#	copy_VarAtts(varin,varout)
#	varout@coordinates = "XLONG XLAT"
#	varout@stagger = " "
#

    if ( unstagDim == "Y") | (unstagDim == "V" ):
	dimV = dims[nd-2]
	if ( nd == 5 ):
	    #varout = 0.5*(varin[:,:,:,:dimV-2,:]+varin[:,:,:,1:dimV-1,:]) #python indexing appears to work differently than ncl, 1 less
	    varout = 0.5*(varin[:,:,:,:dimV-1,:]+varin[:,:,:,1:dimV-0,:])

	if ( nd == 4 ):
	    #varout = 0.5*(varin[:,:,:dimV-2,:]+varin[:,:,1:dimV-1,:])
	    varout = 0.5*(varin[:,:,:dimV-1,:]+varin[:,:,1:dimV-0,:])

	if ( nd == 3 ):
	    #varout = 0.5*(varin[:,:dimV-2,:]+varin[:,1:dimV-1,:])
	    varout = 0.5*(varin[:,:dimV-1,:]+varin[:,1:dimV-0,:])

	if ( nd == 2 ):
	    #varout = 0.5*(varin[:dimV-2,:]+varin[1:dimV-1,:])
	    varout = 0.5*(varin[:dimV-1,:]+varin[1:dimV-0,:])

#	do i = 0,nd-1
#	    varout!i = varin!i
#	end do
#	i = nd-2
#	varout!i = "south_north"
#	copy_VarAtts(varin,varout)
#	varout@coordinates = "XLONG XLAT"
#	varout@stagger = " "


    if ( unstagDim == "Z" ):
	dimW = dims[nd-3]
	if ( nd == 5 ):
	    #varout = 0.5*(varin[:,:,0:dimW-2,:,:]+varin[:,:,1:dimW-1,:,:])
	    varout = 0.5*(varin[:,:,0:dimW-1,:,:]+varin[:,:,1:dimW-0,:,:])

	if ( nd == 4 ):
	    #varout = 0.5*(varin[:,0:dimW-2,:,:]+varin[:,1:dimW-1,:,:])
	    varout = 0.5*(varin[:,0:dimW-1,:,:]+varin[:,1:dimW-0,:,:])

	if ( nd == 3 ):
	    #varout = 0.5*(varin[0:dimW-2,:,:]+varin[1:dimW-1,:,:])
	    varout = 0.5*(varin[0:dimW-1,:,:]+varin[1:dimW-0,:,:])

#	do i = 0,nd-1
#	    varout!i = varin!i
#	end do
#	i = nd-3
#	varout!i = "bottom_top" 
#	copy_VarAtts(varin,varout)
#	varout@coordinates = "XLONG XLAT"
#	varout@stagger = " "

    if ( unstagDim == "M" ):
	print 'staggared dim is "M", no support for this yet'
	varout = varin





    for ele in ["X","U","Y","V","Z", "M"]:
	if (unstagDim == ele):
	    return varout

#    if( any( unstagDim == ["X","U","Y","V","Z"] ) ):
#	return varout  
#    else:
#	print "Warning: Input field was not unstaggered" 
#    return varin 



