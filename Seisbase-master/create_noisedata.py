#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:46:28 2018
convert the miniseed file from IRIS to input formats required by ancc and 
python code
Jan. 20, 2018, use new function merge_seed to create monthly seed 
@author: yunfeng
"""
from obspy import read
import glob
import os
from database import Database,Station,Seed
import pickle 
# load the database object
file_db = open('noise_database.obj', 'r') 
db = pickle.load(file_db)
from par import Parfile
parfile=Parfile()
# define two station that need to be cross-correlated with
net1='IU'
sta1='BRGNA'
#net2='IU'
#sta2 = 'DWPF'
db_sta1 = db.select(network=net1,station=sta1)
#db_sta2 = db.select(network=net2,station=sta2)
output_name = db_sta1.networks[0].code+'.'+db_sta1.networks[0].stations[0].code + '.BHZ.mseed'
for seed1 in db_sta1.networks[0].stations[0].seeds:
    year='{:04d}'.format(seed1.time.year)
    month='{:02d}'.format(seed1.time.month)    
    yyyymm = int(year+month)
    miniseed_directory = parfile.miniseed_directory+year+'-'+month
    # create directory if not exits
    if not os.path.exists(miniseed_directory):
        os.makedirs(miniseed_directory)
    # merge seed with an existing seed in the target directory
    seed1.merge_seed(target_directory=miniseed_directory,target_seed=output_name)
    
