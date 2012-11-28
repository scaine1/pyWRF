#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import numpy as n
import pylab as pl
import pyWRF
import skewt
import dircache
import os

MP='ferrier'
LEVEL=5

#input_directory = path to your wrf files (needs / at the end)
input_directory = '/nfs/monash/home/louwil/WRFV3/031303/wrfout_d03/'

wrf_files=dircache.listdir(input_directory)


for wfile in wrf_files:
    if wfile[0:6] != 'wrfout':
	continue
    
    print wfile
    wfilename=wfile[11:len(wfile)]
    wrf_file=pyWRF.wrf_file(input_directory+wfile)

    
    
    #plot a cappi at model level 10, technically this is not a cappi
    #wrf_file.plot_cappi(level=10,title='model level 10')i


    #return the data into an array so we can work with it
    dbz=wrf_file.calculate_dbz_ferrier()
    
    height_levels=n.arange(500.,20000.,500)   #height of the levels you want in meters
#    print height_levels
    Z=wrf_file.get_var('Z')


    dbz_interp=wrf_file.interp_to_height(height_levels,Z,dbz)

    wrf_file.plot_directory=input_directory+'dbz'

    direxist=os.path.isdir(wrf_file.plot_directory)
    if (direxist == False):
	os.system('mkdir -p '+wrf_file.plot_directory)


    wrf_file.plot_cappi(level=LEVEL,MPscheme=MP,dbz_array=dbz_interp,title=wfilename+' '+str(height_levels[LEVEL]/1000)+'km',imagename=str(height_levels[LEVEL]/1000)+',interp')
    
    

    print 'finished plotting time period '

