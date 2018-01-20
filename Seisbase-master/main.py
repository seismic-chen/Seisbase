#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 23:54:50 2017
scan the CRANE dataset and create the database class
Jan. 19, 2018, scan other networks e.g., RV and TD
@author: yunfeng
"""
import os
import fnmatch
import re
from obspy import UTCDateTime
from database import Database,Network,Station,Seed
import pickle
from par import Parfile

parfile=Parfile()

db = Database(networks=[],code='noise_database')
database_directory=parfile.database_directory
pattern='FullSeed_*.seed'
network_all=next(os.walk(database_directory))[1]
for network_code in network_all:
    if len(network_code) != 2:
        continue
    network_directory=os.path.join(database_directory,network_code)
    net=Network(stations=[],code=network_code,path=network_directory)
    station_all=next(os.walk(network_directory))[1]
    for station_code in station_all:
#        if len(station_code) != 3:
        # change selection criteria, check if all upper case
        if not station_code.isupper():
            continue
        station_directory=os.path.join(database_directory,network_code,station_code)
        # create a station object
        # Note that to initialize a new station class, need to clear seed list
        sta=Station(seeds=[],code=station_code,path=station_directory)
        for root, dirs, files in os.walk(station_directory):
            print root
            for filename in fnmatch.filter(files,pattern):
                path=os.path.join(root,filename)
                datestr=re.split('_|\.',filename)[2]
                timestr=re.split('_|\.',filename)[3]
                se=Seed(code=filename,path=path,time=UTCDateTime(datestr+timestr))
                sta.add(se)
        # add new station to network
        net.add(sta)
    # add new network to database
    db.add(net)
#db.stations[0].get_statics()
t0=UTCDateTime(2007,1,9,1,0,0,0)
t1=t0+100
se = db.networks[0].stations[0].seeds[-2]
st = se.read_seed()
# save the object
object_db = db
file_db = open('noise_database.obj', 'w') 
pickle.dump(object_db, file_db)
file_db.close()
#st.plot()

#se1 = Seed()
#se1.code = 'FullSeed_JOF_20070109_000000.seed'
#se2 = Seed()
#se2.code = 'FullSeed_JOF_20070110_000000.seed'
#st.seed.append(se1)
#st.seed.append(se2)
# output the miniseed file according to requirment by noise tomography structure