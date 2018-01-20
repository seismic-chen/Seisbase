#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 15:44:02 2017

@author: yunfeng
create the station xml for each CRANE station
"""

from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy import UTCDateTime
from obspy import read
from obspy.io.xseed import Parser
from obspy import read_inventory
import matplotlib.pyplot as plt
output_directory = '/Users/yunfeng/20_30/research/python_codes/seismic-noise-tomography-master/StationXML'
inv = read_inventory("CRANE.xml")
for station in inv.networks[0].stations:
    code=station.code
    invtmp=inv.select(station=code)
    filename=code+".xml"
    filepath = os.path.join(output_directory, filename)
    invtmp.write(filepath,format="STATIONXML")
    