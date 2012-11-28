#!/usr/bin/python
import numpy as n
import pylab as pl
import pyWRF
import skewt
import dircache
import os

#input_directory = path to your wrf files (needs / at the end)
input_directory = '/data/cdgs1/scaine/bc_surface_test/spinup/d01/'

wrf_files=dircache.listdir(input_directory)


for wfile in wrf_files:
    if wfile[0:6] != 'wrfout':
	continue

    wrf_file=pyWRF.wrf_file(input_directory+wfile)
    times=wrf_file.get_var('Times')
    u_wind=wrf_file.get_var('U')
    xlat  =wrf_file.get_var('XLAT')
    xlon  =wrf_file.get_var('XLONG')
    melb_airport_lat=-37.67
    melb_airport_lon=144.83


    #note this will fail if your domain does not include melbourne
    #i,j = wrf_file.get_ij_lat_long(xlat,xlon,melb_airport_lat,melb_airport_lon)

    darwin_airport_lat=-12.4078
    darwin_airport_lon=130.0876
    #i,j = wrf_file.get_ij_lat_long(xlat,xlon,darwin_airport_lat,darwin_airport_lon)
    i,j= 0,0

    print i,j
    wrf_file.plot_directory='plots'
    direxist=os.path.isdir(wrf_file.plot_directory)
    if (direxist == False):
	os.system('mkdir -p '+wrf_file.plot_directory)

    for time_count in range(len(times)):
	time=times[time_count]
	wrf_file.plot_skewt(i,j,time+'_skewt.png','Darwin Airport')
        wrf_file.plot_vapor(imagename=time+'_vapor.png',title=time,timestep=time_count)
	print 'finished plotting time period ' + time
