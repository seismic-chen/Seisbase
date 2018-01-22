#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 09:24:07 2017
Database class, contains Database, Station and Seed objects
Jan. 18, 2018, add network class such that the database could read in data from
other network
Jan. 20, 2018, fix the error in select function
add exception when reading miniseed file
@author: yunfeng
"""

import matplotlib.pyplot as plt
from matplotlib.dates import date2num,datetime,DateFormatter
import numpy as np
from obspy import read
import os
import copy
import fnmatch
from par import Parfile
import shutil
parfile=Parfile()

class Database(object):
    """ Database class contains Network class """
    
    def __init__(self, networks=list(),code=''):
        """ Create a new networks class """
        self.networks = networks
    def add(self,networkfile):
        """ add new network class """
        self.networks.append(networkfile)
        
    def select(self,network=None,station=None):
        """ Select network (station) that matches network (station) code, modified after 
        obspy select function"""
        networks = []
        for net in self.networks:
        # skip if any given criterion is not matched
            if network is not None:
                if not fnmatch.fnmatch(net.code.upper(),
                                       network.upper()):
#                    print net.code+'skip'
                    continue
            has_stations = bool(net.stations)
            if not has_stations:
                continue
            
            net_=net.select(station=station)
            # check if has station
            has_stations = bool(net_.stations)
            if not has_stations:
                continue
            networks.append(net_)
        database = copy.copy(self)
        database.networks = networks
        return database
    
class Network(object):
    """ Network class contains the Station class """
    def __init__(self,stations=list(),code='',path='',start_date='',end_date='',
                 total_number_of_stations=None):
        self.stations = stations
        self.code = code
        self.path = path
        self.start_date=str(start_date)
        self.total_number_of_stations=total_number_of_stations
        self.end_date=str(end_date)
        
    def select(self,station=None):
        stations = []
        # skip if any given criterion is not matched
        for sta in self.stations:
            if station is not None:
                if not fnmatch.fnmatch(sta.code.upper(),
                                       station.upper()):    
                    continue
            stations.append(sta)
        networks = copy.copy(self)
        networks.stations = stations
        return networks
        
    def add(self,stationfile):
        """ add new Seed class """
        self.stations.append(stationfile)
        
    def __str__(self):
        """ Modified after obspy __str__ function """
        contents = self.get_contents()     
        ret = ("Network {network_name}\n"
               "\Network Code: {station_code}\n"
               "\tStation Count: A total of {total} stations\n"
               "\tStart-End: {start_date} - {end_date}\n")
        ret = ret.format(
            network_name=contents["networks"][0],
            station_code=self.code,
            total=self.total_number_of_stations,
            start_date=str(self.start_date) if self.start_date else "",
            end_date=str(self.end_date) if self.end_date else "")
        return ret
    
    def _repr_pretty_(self, p, cycle):
        p.text(str(self))
        
    def get_contents(self):
        """ Modified after obspy get_contents function. """
        desc = "%s" % (self.code)
        content_dict = {"networks": [desc], "stations": []}
        self.get_total_number_of_stations()
        
        for station in self.stations:
            content_dict["stations"].append(
                "%s.%s.%s" % (self.code, station.code, station.start_date))
        return content_dict   
     
    def get_total_number_of_stations(self):
        """ count the total number of seeds file """
        self.total_number_of_stations = len(self.stations)
        return self.total_number_of_stations    
    
class Station(object):
    """ Station class contains the Seed class """
    def __init__(self,seeds=list(),code='',path='',start_date='',end_date='',
                 latitude=None,longitude=None,elevation=None,total_number_of_seeds=None,
                 dateless=None):
        self.seeds = seeds
        self.code = code
        self.path = path
        self.start_date=str(start_date)
        self.total_number_of_seeds=total_number_of_seeds
        self.end_date=str(end_date)
        self.latitude=latitude
        self.longitude=longitude
        self.elevation=elevation
        self.dateless=dateless
#    def __repr__(self):
#    def __str__(self):

    def get_statics(self):
        import matplotlib.cm as cm
        """ Return data statistics for a given station """
        # prepare date string list
        time_string_list = list()
        yday_list = list()
        year_list = list()
        for n in range(0,len(self.seeds)):
            tmp_time_string = self.seeds[n].time
            time_string_list.append(tmp_time_string.datetime)
            year_list.append(tmp_time_string.timetuple().tm_year)
            yday_list.append(tmp_time_string.timetuple().tm_yday)
        #x = np.array(date2num(time_string_list))
        #y = np.ones(40)
        #plt.plot_date(x,y)
        #t0 = datetime.date(2006,01,01)
        #t1 = datetime.date.today()
        #fig, ax = plt.subplots(1)
        #fig.autofmt_xdate()
        #plt.plot_date(x, y)
        #xfmt = DateFormatter('%d-%m-%y')
        #ax.xaxis.set_major_formatter(xfmt)
        #ax.set_xlim([t0,t1])
        #ax.set_ylim([0, 2])
        #plt.show()
        # a chart similar to the one on Github
        # convert UTCDateTime to datetime string
        t_min = min(time_string_list)
        t_max = max(time_string_list)
        # find the Sunday before t_min and Saturday after t_max
        t_start = t_min - datetime.timedelta(days=t_min.weekday())
        t_end = t_max + datetime.timedelta(days=6-t_max.weekday())
        time_range = (t_end-t_start).days+1
        # recreate time axis
        dt = t_start
        step = datetime.timedelta(days=1)
        time_axis = []
        while dt <= t_end:
            time_axis.append(dt.strftime('%Y-%m-%d'))
            dt += step
        # compute the time difference in days
        time_diff = []
        time_diff = [(x - t_start).days for x in time_string_list]
        # convert list to numpy array
        keep_index = np.asarray(time_diff)
        # create a time vector
        C = np.zeros(time_range)
        C[keep_index] = 1
        T = np.asarray(time_axis)
        # reshape C to a matrix such that each column represents a week
        A = C.reshape((7,len(C)/7),order='F')
        B = T.reshape((7,len(C)/7),order='F')
        # get the shape of A
        [M,N] = np.shape(A)
        plt.figure(figsize=(18,10))
        plt.imshow(A,vmin=0,vmax=3,cmap=cm.Greens,interpolation='none',aspect='equal')
        plt.grid(True)
        ax = plt.gca()
        
#        labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
#                  'Sep', 'Oct', 'Nov', 'Dec']
        ax.set_xticks(np.arange(-.5, N, 1))
        ax.set_yticks(np.arange(-.5, M, 1))
        # ax.set_xticklabels(labels)
        ylabels = ('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
        ax.set_yticklabels(ylabels)
        
        x = np.arange(-.5, N, 1)
        xlabels = B[0,:]
        # plot label every 4 ticks
        indices = range(0,len(xlabels),4)
        mask = np.ones(xlabels.shape,dtype=bool)
        mask[indices] = False
        xlabels[mask] = ''
        # You can specify a rotation for the tick labels in degrees or with keywords.
        plt.xticks(x, xlabels, rotation='vertical')
        ax.xaxis.tick_top()
        
    def get_inventory(self):
        from obspy import read_inventory
        inventory_path = parfile.inventory_path
        inv = read_inventory(inventory_path,format='STATIONXML')
        network = inv.networks[0].select(station=self.code)
        self.dateless = network
        self.latitude = network[0][0].latitude
        self.longitude = network[0][0].longitude
        self.elevation = network[0][0].elevation
        return self
    
    def select(self,time=None):
        """ Select seed that matches the given time, modified after 
        obspy select function"""
        seeds = list()
        for seed in self.seeds:
        # skip if any given criterion is not matched
            if seed is not None:
                if not seed.time == time:
                    continue
            seeds.append(seed)
        station = copy.copy(self)
        station.seeds = seeds
        return station
        
    def add(self,seedfile):
        """ add new Seed class """
        self.seeds.append(seedfile)
        
    def __str__(self):
        """ Modified after obspy __str__ function """
        contents = self.get_contents()
        self.get_total_number_of_seeds()
        
        ret = ("Station {station_name}\n"
               "\tStation Code: {station_code}\n"
               "\tSeed Count: A total of {total} seeds\n"
               "\tStart-End: {start_date} - {end_date}\n"
               "\tLatitude: {lat:.2f}, Longitude: {lng:.2f}, "
               "Elevation: {elevation:.1f} m\n"
               "\tDateless: {dateless}\n")
        ret = ret.format(
            station_name=contents["stations"][0],
            station_code=self.code,
            total=self.total_number_of_seeds,
            start_date=str(self.start_date) if self.start_date else "",
            end_date=str(self.end_date) if self.end_date else "",
            lat=self.latitude if self.latitude else -999, 
            lng=self.longitude if self.longitude else -999,
            elevation=self.elevation if self.elevation else -999,
            dateless=self.dateless if self.dateless else '')
        return ret
    
    def _repr_pretty_(self, p, cycle):
        p.text(str(self))
        
    def get_contents(self):
        """ Modified after obspy get_contents function. """
        desc = "%s" % (self.code)
        content_dict = {"stations": [desc], "seeds": []}

        for seed in self.seeds:
            content_dict["seeds"].append(
                "%s.%s.%s" % (self.code, seed.code, seed.time))
        return content_dict
    
    def get_total_number_of_seeds(self):
        """ count the total number of seeds file """
        self.total_number_of_seeds = len(self.seeds)
        return self.total_number_of_seeds
        
class Seed(object):
    """ Seed class contains the FullSeed file information """
    def __init__(self,code='',path='',time=''):
        self.code = code
        self.path = path
        self.time = time
        
    def read_seed(self,t0=None,t1=None):
        """ Read seed file using obspy read module
        Obspy seems only support miniseed format, so a walkaround is to convert
        Fullseed file to miniseed using rdseed first, then read in with obspy
        """
        import subprocess
        # read the miniseed file
        try:
            if t0 and t1:
                st = read(self.path,starttime=t0, endtime=t1)
            else:
                st = read(self.path)
        except:
            cwd = os.getcwd()
            seed_directory = os.path.dirname(self.path)
            file_name = os.path.basename(self.path)        
            os.chdir(seed_directory)
            rdseed_cmd = parfile.rdseed_path    
            # if miniseed file exists, remove the file first
            if os.path.exists('mini.seed'):
                os.unlink('mini.seed')
            subprocess.call([rdseed_cmd, '-d', '-o', '4', '-f', file_name])
            if t0 and t1:
                st = read('mini.seed',starttime=t0, endtime=t1)
            else:
                st = read('mini.seed')
            # delete miniseed
            os.unlink('mini.seed')
            os.chdir(cwd)
        return st
    
    def move_seed(self,output_directory=None,output_name=None,
                            channel=None):
        """ Move (make a copy of) a daily seed/miniseed files to the target 
        direcotry        
        """
        import subprocess

        filename=os.path.basename(self.path)
        default_filename = os.path.join(output_directory, filename)
        if output_name is not None:
            default_filename = os.path.join(output_directory, output_name)
        self.miniseed_direcotry=output_directory
        shutil.copyfile(self.path,default_filename)
        return self.miniseed_direcotry
    
    def merge_seed(self,target_directory=None,target_seed=None):
        """ Merge a daily seed/miniseed files with an target seed
        """
        import subprocess
        if not os.path.exists(self.path):  
            print "seed file doesn't exists"  
            raise IOError  
        filein=self.path
        fileout=os.path.join(target_directory, target_seed)
        # check if output file exists
        if os.path.isfile(fileout):
            mode='ab'
        else:
            mode='wb'           
        with open(fileout,mode) as output: 
            with open(filein, 'rb') as input:  
                data=input.read()  
                output.write(data)  
    
    def convert_to_miniseed(self,output_directory=None,output_name=None,
                            channel=None):
        """ Convert the Fullseed file to miniseed
        could append argument to extract certain station or component
          -C arg retrieve the comments where 'arg' is either STN or CHN"""
        import subprocess   
        cwd = os.getcwd()
        seed_directory = os.path.dirname(self.path)
        file_name = os.path.basename(self.path)
        os.chdir(seed_directory)
        if output_directory is None:
            output_directory = cwd
        rdseed_cmd = parfile.rdseed_path
        subprocess.call([rdseed_cmd, '-d', '-o', '4', '-f', file_name, '-q',output_directory])
        # check if miniseed exists
        default_filename = os.path.join(output_directory, "mini.seed")
        if not os.path.isfile(default_filename):
            print "Failed to read seed"
            return None 
        # if miniseed file exists
        if output_name is not None:
            new_filename = os.path.join(output_directory, output_name)
            os.rename(default_filename, new_filename)
            self.miniseed_path = new_filename
        else:
            self.miniseed_path = default_filename
#        if t0 and t1:
#            st = read('mini.seed',starttime=t0, endtime=t1)
#        else:
#            st = read('mini.seed')
        os.chdir(cwd)
        return self.miniseed_path
    
    def __str__(self):
        """ Modified after obspy __str__ function """        
        ret = ("Seed {seed_code}\n"
               "\tSeed Code: {seed_code}\n"
               "\tDate: {date}\n"
               "\tPath: {path}\n")
        ret = ret.format(
            seed_code=self.code,
            date=str(self.time) if self.time else "",
            path=self.path)
        return ret
    
    def _repr_pretty_(self, p, cycle):
        p.text(str(self))
#    def request_instrument_response
#    def distance_from_orogin(self):
#        """ computer my distance from the origin """
#        return ((self.x ** 2) + (self.y ** 2)) ** 0.5
#    def __add__(self, other):
#        new = copy.deepcopy(self)
#        new += other
#        return new
