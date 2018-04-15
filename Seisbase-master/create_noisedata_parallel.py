#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:46:28 2018
convert the miniseed file from IRIS to input formats required by ancc and 
python code
Jan. 20, 2018, use new function merge_seed to create monthly seed
Mar. 23, 2018, set the time range for data conversion; convert FullSeed to 
miniseed for CRANE station
Mar. 24, 2018, make the code parallel using multiprocessing module
Apr. 12, 2018, code has bug when extracting Z component using parallel
@author: yunfeng
"""
from obspy import read,UTCDateTime
import glob
import os
from database import Database,Station,Seed
import pickle 
import multiprocessing as mp
# load the database object
file_db = open('noise_database.obj', 'r') 
db = pickle.load(file_db)
from par import Parfile
parfile=Parfile()
# define two station that need to be cross-correlated with
#net1='RV'
#sta1='DWPF'
#net2='IU'
#sta2 = 'DWPF'
start_time=UTCDateTime(2016,1,8)
end_time=UTCDateTime(2016,3,31)
#db_sta1 = db.select(network=net1,station=sta1)
#db_net1 = db.select(network=net1)
nseed=int(end_time-start_time)/(3600*24)+1
#station=db_net1.networks[0].stations[0]
#seed1=station.seeds[0]
#miniseed_name = station.network_code+'.'+station.code + '.BHZ.mseed'
#year='{:04d}'.format(seed1.time.year)
#month='{:02d}'.format(seed1.time.month)    
#miniseed_directory = parfile.miniseed_directory+year+'-'+month+'/'+station.code  
#miniseed_path=os.path.join(miniseed_directory,miniseed_name)
#seed1.convert_to_miniseed(output_directory=miniseed_directory,output_name=miniseed_name,channel='*Z',if_merge=True)

def seed_convertor(station,start_time,end_time):
    """ function converts the seed file"""
    proc = os.getpid()
    print('Converting Station {0} by process id: {1}'.format(
        station.code, proc))
    miniseed_name = station.network_code+'.'+station.code + '.BHZ.mseed'
    k=0
    for seed1 in station.seeds:
        time=seed1.time
        if time>=start_time and time<=end_time:
            k=k+1
            print('Station {0:s}: Merging {1:d} out of {2:d} seeds'.format(station.code,k,nseed))
            year='{:04d}'.format(seed1.time.year)
            month='{:02d}'.format(seed1.time.month)    
            miniseed_directory = parfile.miniseed_directory+year+'-'+month+'/'+station.code            
            # create directory if not exits
            if not os.path.exists(miniseed_directory):
                os.makedirs(miniseed_directory)
            miniseed_path=os.path.join(miniseed_directory,miniseed_name)
            # convert Fullseed to miniseed for CR network
            if station.network_code == 'CR' or station.network_code == 'CN' :
                seed1.convert_to_miniseed(output_directory=miniseed_directory,output_name=miniseed_name,if_merge=True,channel='*Z')
            else:
                # merge seed with an existing seed in the target directory
                seed1.merge_seed(target_directory=miniseed_directory,target_seed=miniseed_name,channel='*Z')
    return None
for network in db.networks[2:4]:
    for station in network.stations[1:]:
        seed_convertor(station,start_time, end_time)

#for station in db_net1.networks[0].stations[0:2]:
#    procs = []
#    proc = mp.Process(target=seed_convertor, args=(station,start_time,end_time,))
#    procs.append(proc)
#    proc.start()
# 
#for proc in procs:
#    proc.join()

#db_sta2 = db.select(network=net2,station=sta2)
#miniseed_name = db_sta1.networks[0].code+'.'+db_sta1.networks[0].stations[0].code + '.BHZ.mseed'
#k=0
#nseed=len(db_sta1.networks[0].stations[0].seeds)
#for seed1 in db_sta1.networks[0].stations[0].seeds:
#    k=k+1
#    print('Merging {0:d} out of {1:d} seeds'.format(k,nseed))
#    year='{:04d}'.format(seed1.time.year)
#    month='{:02d}'.format(seed1.time.month)    
#    yyyymm = int(year+month)
#    miniseed_directory = parfile.miniseed_directory+year+'-'+month
#    # create directory if not exits
#    if not os.path.exists(miniseed_directory):
#        os.makedirs(miniseed_directory)
#    miniseed_path=os.path.join(miniseed_directory,miniseed_name)
#    # merge seed with an existing seed in the target directory
#    seed1.merge_seed(target_directory=miniseed_directory,target_seed=miniseed_name)