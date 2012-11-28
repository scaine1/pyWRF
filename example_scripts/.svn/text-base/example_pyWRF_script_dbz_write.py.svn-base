#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import numpy as n
import pylab as pl
import pyWRF
import skewt
import dircache
import os
import cPickle as pickle

#input_directory = path to your wrf files (needs / at the end)
input_directory = '/nfs/monash/home/louwil/WRFV3/022403_Lin/wrfout_d03/'

wrf_files=dircache.listdir(input_directory)


for wfile in wrf_files:
    if wfile[0:6] != 'wrfout':
	continue

    wrf_file=pyWRF.wrf_file(input_directory+wfile)
    times=wrf_file.get_var('Times')


    #plot a cappi at model level 10, technically this is not a cappi
    #wrf_file.plot_cappi(level=10,title='model level 10')


    #return the data into an array so we can work with it
    dbz=wrf_file.calculate_dbz_lin()
    height_levels=n.arange(500.,20000.,500)   #height of the levels you want in meters
#    print height_levels
    Z=wrf_file.get_var('Z')


    dbz_interp=wrf_file.interp_to_height(height_levels,Z,dbz)

    wrf_file.plot_directory= input_directory+'dbz_data'
    direxist=os.path.isdir(wrf_file.plot_directory)
    if (direxist == False):
	os.system('mkdir -p '+wrf_file.plot_directory)

    # Example: 4 dimensions ['T','BT','SN','WE']
    wrf_file.write_netcdf_file(input_variable=dbz_interp[:,:,:,:],var_name='dBZ_interp',directory=wrf_file.plot_directory,filename=times[0]+'_ex1_dbz',var_dim=('T','BT','SN','WE'))  
    wrf_file.write_netcdf_file(input_variable=times,var_name='Times',directory=wrf_file.plot_directory,filename=times[0]+'_ex1_dbz',var_dim=('T'))


#   for the_time in range(len(times)):
 	
	# Example: 3 dimensions ['BT','SN','WE']
#   	wrf_file.write_netcdf_file(input_variable=dbz_interp[the_time,:,:,:],var_name='dBZ_interp',directory=wrf_file.plot_directory,filename=times[the_time]+'_ex2_dbz',var_dim=('BT','SN','WE'))  
 	
	# Example: 3 dimensions ['T','SN','WE']   	
#	wrf_file.write_netcdf_file(input_variable=dbz_interp[:,0,:,:],var_name='dBZ_interp',directory=wrf_file.plot_directory,filename=times[the_time]+'_ex3_dbz',var_dim=('T','SN','WE'))  
    break
#    wrf_file.plot_cappi(level=10,dbz_array=dbz_interp,title=height_levels[level]/10)+'km',imagename='interp')
       
    

#    print 'finished plotting time period '

