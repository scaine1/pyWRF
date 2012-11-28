#!/usr/bin/python
import numpy as n
import pylab as pl
import pyWRF
import skewt
import dircache
import os

#input_directory = path to your wrf files (needs / at the end)
input_directory = '/data/cdgs1/scaine/bc_surface_test3/tracers_added/d01/'

wrf_files=dircache.listdir(input_directory)

first_time=True
for wfile in wrf_files:
    if wfile[0:6] != 'wrfout':
	continue
    if (first_time== False):
	continue

    wrf_file=pyWRF.wrf_file(input_directory+wfile)
    times=wrf_file.get_var('Times')
    u_wind=wrf_file.get_var('U')
    v_wind=wrf_file.get_var('V')

    wrf_file.write_netcdf_file(u_wind,'u_wind',filename='blah',directory='netcdf_out')

    wrf_file.write_netcdf_file(v_wind,'v_wind',filename='blah',directory='netcdf_out')

    first_time=False
