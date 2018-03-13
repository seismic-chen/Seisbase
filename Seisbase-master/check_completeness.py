#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 10:54:06 2018
Mar. 4, 2018, Y.C., Code to check the completeness of the downloaded ambient 
noise data
@author: yunfeng
"""
import os.path
from obspy.clients.fdsn.client import Client
from obspy import UTCDateTime
from obspy import read, read_inventory
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from database import Database,Station,Seed
import pickle 

inv_crane = read_inventory('CRANE.xml')
#inv_crane = read_inventory('RV.xml')

starttime=UTCDateTime("2014-01-01")
endtime=UTCDateTime("2017-12-31")
client = Client("IRIS")
inventory = client.get_stations(starttime=starttime,endtime=endtime,
                                minlatitude=47.5,maxlatitude=60.5,
                                minlongitude=-122,maxlongitude=-108)
inv_raven = inventory.select(network='RV')
inv_td = inventory.select(network='TD')
inv_cnsn = inventory.select(network='CN')

fig, ax = plt.subplots(figsize=(12, 8))
m=Basemap(width=1000000,height=1500000,
            resolution='l',projection='eqdc',\
            lat_1=50.,lat_2=60,lat_0=55,lon_0=-115.)
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcoastlines()
m.drawcountries()
m.drawstates()
# Draw a lon/lat grid (20 lines for an interval of one degree)
m.drawparallels(np.linspace(48, 60, 7), labels=[1, 1, 0, 0], fmt="%.2f",
                dashes=[2, 2])
m.drawmeridians(np.linspace(-120, -110, 2), labels=[0, 0, 1, 1], fmt="%.2f",
                dashes=[2, 2])
fig.bmap = m
colors = {'CN': 'blue', 'TD': 'green', 'RV': 'red', 'CR': 'purple'}
inv_raven.plot(fig=fig, show=False,label=False,color_per_network=colors)
inv_td.plot(fig=fig, show=False,label=False,color_per_network=colors)
inv_cnsn.plot(fig=fig, show=False,label=False,color_per_network=colors)
inv_crane.plot(fig=fig, show=False,label=False,color_per_network=colors)
plt.show()

# load the database object
from database import Database,Station,Seed
import pickle 
file_db = open('noise_database.obj', 'r') 
db = pickle.load(file_db)
#stations = db.networks[0].stations
#station=stations[0]

networks = db.networks
start_date=UTCDateTime(2015,1,1)
end_date=UTCDateTime(2016,12,31)
for network in networks[0:1]:
    stations= network.stations
    for station in stations:
        print station.code
        # check if output exist, if so, skip the current station
        filename=station.code+'_missing_data.txt'
        if os.path.isfile(filename):
            continue
        station.check_completeness(start_date=start_date, end_date=end_date, output_to_file=False)
