#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 13:29:23 2017
@author: yunfeng
# python code to create the input data
# Nov 13, 2016, Yunfeng Chen, Global Seismology Group, University of Alberta
# Nov. 14, 2017, Y.C., GSG, UofA.
# Use database class and change the code to reflect more general use

The required format is:
In the current version of the program, miniseed files MUST be
# organized inside their directory as:
# <year>-<month>/<network>.<station>.<channel>.mseed, e.g.:
# 1988-10/BL.JFOB.BHZ.mseed
# So, there is one sub-directory per month, and inside it, one miniseed
# file per month and per station.
"""

from obspy import read
import glob
import os
from database import Database,Station,Seed
import pickle 
# load the database object
file_db = open('CR_database.obj', 'r') 
db = pickle.load(file_db)

# define two station that need to be cross-correlated with
sta1 = 'BRU'
#sta2 = 'CLA'
db_sta1 = db.select(station=sta1)
#db_sta2 = db.select(station=sta2)
output_name = 'CR.'+ db_sta1.stations[0].code + '.BHZ.mseed'
yyyymm_old = 0
for seed1 in db_sta1.stations[0].seeds:
    year='{:04d}'.format(seed1.time.year)
    month='{:02d}'.format(seed1.time.month)    
    yyyymm_new = int(year+month)
    output_directory = '/Users/yunfeng/20_30/research/python_codes/'+\
    'seismic-noise-tomography-master/miniseed_dir/'+year+'-'+month
    # create directory if not exits
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # write the file to mini format
    if yyyymm_old == 0:
        yyyymm_old = yyyymm_new
        
    if yyyymm_old < yyyymm_new:
        yyyymm_old = yyyymm_new
        # check if mini.seed file exists, if exists rename
        if os.path.exists(default_filename):
            os.rename(default_filename, new_filename)
    else:
        seed1.convert_to_miniseed(output_directory=output_directory,output_name='mini.seed')
        default_filename = os.path.join(output_directory, "mini.seed")
        new_filename = os.path.join(output_directory, output_name)
        
# rename the last mini.seed file
os.rename(default_filename, new_filename)