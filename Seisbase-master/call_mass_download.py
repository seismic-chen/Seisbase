#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 10:42:23 2018
excute mass_download if mseed does not exist in the database. 
@author: yunfeng
"""
# load the database object
from obspy import UTCDateTime
from obspy.clients.fdsn.mass_downloader import RectangularDomain, \
    Restrictions, MassDownloader
    
from database import Database,Station,Seed
import pickle
from datetime import timedelta
file_db = open('noise_database.obj', 'r') 
db = pickle.load(file_db)
network=db.select(network='RV').networks[0]

for station in network.stations:
    start_date=UTCDateTime(2015,1,1)
    end_date=UTCDateTime(2016,12,31)
    station.check_completeness(start_date=start_date, end_date=end_date, output_to_file=True)
    missing_date=station.missing_date
    station_code=station.code
    network_code=station.network_code
    # Rectangular domain containing parts of southern Germany.
    #domain = RectangularDomain(minlatitude=51, maxlatitude=56,
    #                           minlongitude=-120.0, maxlongitude=-112.0)
    for date in missing_date:
        tstart=date
        tend=tstart+timedelta(days=1)
        # Alberta
        domain = RectangularDomain(minlatitude=51, maxlatitude=56,
                                   minlongitude=-120.0, maxlongitude=-112.0)
        restrictions = Restrictions(
            # Get data for a whole year.
            starttime=tstart,
            endtime=tend,
            # Chunk it to have one file per day.
            chunklength_in_sec=86400,
            # Considering the enormous amount of data associated with continuous
            # requests, you might want to limit the data based on SEED identifiers.
            # If the location code is specified, the location priority list is not
            # used; the same is true for the channel argument and priority list.
            network=network_code, station=station_code, location="*", channel="*",
            # The typical use case for such a data set are noise correlations where
            # gaps are dealt with at a later stage.
            reject_channels_with_gaps=False,
            # Same is true with the minimum length. All data might be useful.
            minimum_length=0.0,
            # Guard against the same station having different names.
            minimum_interstation_distance_in_m=100.0)
    
        # Restrict the number of providers if you know which serve the desired
        # data. If in doubt just don't specify - then all providers will be
        # queried.
        mdl = MassDownloader(providers=["IRIS"])
        mdl.download(domain, restrictions, mseed_storage="waveforms",
                     stationxml_storage="stations")


#def get_mseed_storage(network, station, location, channel, starttime,
#                       endtime):
#     # Returning True means that neither the data nor the StationXML file
#     # will be downloaded.
#     filename="{}"
#     if is_in_db(network, station, location, channel, starttime, endtime):
#         return True
#     # If a string is returned the file will be saved in that location.
#     return os.path.join(ROOT, "%s.%s.%s.%s.mseed" % (network, station,
#                                                      location, channel))
#     
#def is_in_db():
#    filename=
#    if os.path.exists(os.path.dirname(interval.filename)):
#        return True
#    
#    
#    if callable(str_or_fct):
#        path = str_or_fct(network, station, location, channel, starttime,
#                          endtime)
#    elif ("{network}" in str_or_fct) and ("{station}" in str_or_fct) and \
#            ("{location}" in str_or_fct) and ("{channel}" in str_or_fct) and \
#            ("{starttime}" in str_or_fct) and ("{endtime}" in str_or_fct):
#        path = str_or_fct.format(
#            network=network, station=station, location=location,
#            channel=channel, starttime=starttime.strftime(strftime),
#            endtime=endtime.strftime(strftime))
#    else:
#        path = os.path.join(
#            str_or_fct,
#            "{network}.{station}.{location}.{channel}__{s}__{e}.mseed".format(
#                network=network, station=station, location=location,
#                channel=channel, s=starttime.strftime(strftime),
#                e=endtime.strftime(strftime)))
#
#    if path is True:
#        return True
#    elif not isinstance(path, (str, bytes)):
#        raise TypeError("'%s' is not a filepath." % str(path))
#    return path

