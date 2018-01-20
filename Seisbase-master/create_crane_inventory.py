#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 22:10:40 2017

@author: yunfeng
This code reads in existing dateless seed, updates the instrument response 
using new module NRL of obspy and writes to a inventory file contains all 
station information.
"""

# create CRANE inventory file
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.clients.nrl import NRL
from obspy import UTCDateTime
from obspy import read
from obspy.io.xseed import Parser
from obspy import read_inventory
import matplotlib.pyplot as plt

#inv = Inventory(networks=[], 
#                source="Global Seismology Group, University of Alberta")
#
#net = Network(code="CR",stations=[],description="CRANE stations")
#
#sta = Station(code="BR2",
#              latitude=53.3294,
#              longitude=-117.8606,
#              elevation=1030.0,4
#              creation_date=UTCDateTime(year=2013,julday=300),
#              site=Site(name="BR2"))
#
#cha = Channel(code="BHZ",
#              location_code="",
#              latitude=53.3294,
#              longitude=-117.8606,
#              elevation=1030.0,
#              depth=0.0,
#              azimuth=0.0,
#              dip=-90.0,
#              sample_rate=40)
#
#sensor_keys=['Nanometrics', 'Trillium 120P']
#datalogger_keys=['Nanometrics', 'Taurus','40 Vpp (0.4)', 'Low (default)', 'Off', '40']
#nrl = NRL()
#response = nrl.get_response(sensor_keys=sensor_keys,datalogger_keys=datalogger_keys)
#
#cha.response = response
#sta.channels.append(cha)
#net.stations.append(sta)
#inv.networks.append(net)

# load BR2.xml file
sp=Parser("/Users/yunfeng/20_30/study/instrument_list/PDCC/crane.dataless_editted_Nov10_2017")
sp.write_xseed("dataless.seed.CR.xml")

inv = read_inventory("dataless.seed.CR.xml")
network = inv.networks[0]
for i in range(0,len(network)):
    station = network.stations[i]
    # extract response from all channels
#    responses = [x.response for x in channels]
#    indices = [int(i) for i,x in enumerate(responses) if x]
#    channel_templates = [channels[i] for i in indices]  
    channel_templates = [x for x in station.channels if x.response]
    # create response file for BHN and BHZ channels
    for j in range(0,len(channel_templates)):
        channel_BHZ = channel_templates[j].copy()
        channel_BHZ.code = 'BHZ'
        channel_BHZ.azimuth = 0.
        channel_BHZ.dip = -90.0
        
        channel_BHN = channel_templates[j].copy()
        channel_BHN.code = 'BHN'
        channel_BHN.azimuth = 0.
        channel_BHN.dip = 0.
        
        station.channels.append(channel_BHZ)
        station.channels.append(channel_BHN)
    network.stations[i] = station
# combine MCL results
MCL = network.select(station='MCL')
# extract channel object
channel_tmp = [x for sublist in MCL.stations for x in sublist.channels]
#flat_list = [item for sublist in channel_tmp for item in sublist]
#remove all MCL from network object
for n in range(0,len(MCL)):
    network.stations.remove(MCL[n])
## add new MCL station
MCL[0].channels = channel_tmp
network.stations.append(MCL[0])
inv.networks[0] = network
inv.write("CRANE.xml",format="STATIONXML")
# inv1.plot(projection='local')
#cha1 = inv1[0].stations[0].channels[0]
#cha1.response.plot(min_freq=0.01,output='VEL')
#cha.response.plot(min_freq=0.01,output='VEL')
# loop over all stations



#nrl = NRL()
#sensor=nrl.sensors['Nanometrics']['Trillium 120P']
#datalogger=nrl.dataloggers['Nanometrics']['Taurus']['40 Vpp (0.4)']['Low (default)']['Off']['20']
#sensor_keys=['Nanometrics', 'Trillium 120P']
#datalogger_keys=['Nanometrics', 'Taurus','40 Vpp (0.4)', 'Low (default)', 'Off', '20']
#response = nrl.get_response(sensor_keys=sensor_keys,datalogger_keys=datalogger_keys)
#response.plot(0.001,output='VEL')
## read dataless file

## prepare the station file
#sp=Parser("/Users/yunfeng/20_30/study/instrument_list/PDCC/JOF.dataless")
