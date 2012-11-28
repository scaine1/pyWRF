#!/usr/bin/python
import numpy as n
import pylab as pl
import pyWRF
import skewt
import dircache
import os

#input_directory = path to your wrf files (needs / at the end)
input_directory = '/data/cdgs2/datasets/dry_line_case_1/send01/d01/'

wrf_files=dircache.listdir(input_directory)


for wfile in wrf_files:
    if wfile[0:6] != 'wrfout':
	continue

    wrf_file=pyWRF.wrf_file(input_directory+wfile)
    times=wrf_file.get_var('Times')

    wrf_file.plot_directory='plots'
    direxist=os.path.isdir(wrf_file.plot_directory)
    if (direxist == False):
	os.system('mkdir -p '+wrf_file.plot_directory)

    for time_count in range(len(times)):
	time=times[time_count]
	wrf_file.plot_precip(imagename=time+'_precip.png',title=time,timestep=time_count)
	print 'finished plotting time period ' + time
