#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 09:24:07 2017
Database class, contains Database, Station and Seed objects
Jan. 18, 2018, add network class such that the database could read in data from
other network
Jan. 20, 2018, fix the error in select function
add exception when reading miniseed file
Mar. 4, 2018, add function to check the completeness of the database
@author: yunfeng
"""

import matplotlib.pyplot as plt
from matplotlib.dates import date2num,datetime,DateFormatter
import numpy as np
from obspy import read, UTCDateTime
import os
import copy
import fnmatch
from par import Parfile
import shutil
import warnings
import datetime
import csv
from collections import OrderedDict

parfile=Parfile() 

class Database(object):
    """ Database class contains Network class """
    
    def __init__(self, networks=list(),code=''):
        """ Create a new networks class """
        self.networks = networks
    def add(self,networkfile):
        """ add new Network class """
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
    
#   def get_seed_storage(network=None, station=None, location=None, channel=None, 
#                    starttime=None, endtime=None):
#       """ Mar 8. 2018, Y.C., returning True if data exist otherwise return the
#       (excepted) path of the data following the name convension
#       """
       
        
   
class Network(object):
    """ Network class contains the Station class
    Change log: Mar. 4, 2018, Y.C., include network code in station class
    """
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
        """ add new Station class """
        self.stations.append(stationfile)
        self._get_network_code()
        
        
    def __str__(self):
        """ Modified after obspy __str__ function """
        contents = self.get_contents()     
        ret = ("Network {network_name}\n"
               "\tNetwork Code: {station_code}\n"
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
    
    def _get_network_code(self):
        """" get the network code and assign to station class """
        for n in range(0,len(self.stations)):
            self.stations[n].network_code = self.code
        return self
    
    def _create_data_completeness_table(self):
        """Create data completness table following the AGS convension
        e.g., 
        stn	2010 2011 2012 2013   2014 2015 2016 2017		
        BR2	N/A	 N/A  N/A  0.167	 N/A	  N/A  N/A N/A"""
        stations=self.stations
        dicttmp=[x.data_completeness for x in stations]
        filename=self.code+'_data_completeness.csv'
        with open(filename, 'w') as csvfile:
            fieldnames=dicttmp[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(dicttmp)

    
class Station(object):
    """ Station class contains the Seed class 
    Change Log: Mar. 4, 2018, Y.C., add check_completness function
    the data
    """
    def __init__(self,seeds=list(),network_code='',code='',path='',start_date='',end_date='',
                 latitude=None,longitude=None,elevation=None,total_number_of_seeds=None,
                 dataless=None):
        self.seeds = seeds
        self.network_code=network_code
        self.code = code
        self.path = path
        self.start_date=str(start_date)
        self.total_number_of_seeds=total_number_of_seeds
        self.end_date=str(end_date)
        self.latitude=latitude
        self.longitude=longitude
        self.elevation=elevation
        self.dataless=dataless
#    def __repr__(self):
#    def __str__(self):
        
    def check_completeness(self,start_date=None,end_date=None,output_to_file=False):
        """ Return the percentage of data and missing date based on the operation
        period of stations or given start_date and end_date
        """
        # prepare date string list
        time_string_list = list()
        size_string_list = list()
        seed_list = list()
        missing_date = list()
        small_date = list()
        small_size = list()
        for n in range(0,len(self.seeds)):
            tmp_time_string = self.seeds[n].time
            time_string_list.append(tmp_time_string)
            tmp_size_string = float(self.seeds[n].size)/1024./1024. # convert to MB
            size_string_list.append(tmp_size_string)
 
        # obtain the station operation period from stationxml file
        self.get_inventory()
        station = self.dataless
        if start_date is None:
            start_date=station.start_date
        if end_date is None:
            end_date=station.end_date
        self.start_date=start_date
        self.end_date=end_date
        duration=end_date.datetime-start_date.datetime+datetime.timedelta(days=1)
        # create consecutive date list
        numdays=duration.days
        datelist = [start_date + datetime.timedelta(days=x) for x in range(0, numdays)]
                
        if time_string_list:
            if not station.start_date:
                warnings.warn('Station start date is not defined')
                station.start_date = UTCDateTime(1900,1,1,0,0,0)
            else:
                if min(time_string_list) < station.start_date:
                    warnings.warn('Data begins before the first day of the station')
            if not station.end_date:
                warnings.warn('Station end date is not defined')
                station.end_date = UTCDateTime(2599,12,31,23,59,59)
            else:
                if  max(time_string_list) >station.end_date:
                    warnings.warn('Data ends after the last day of the station')
                    
        # keep the seeds in the given time period 
        seed_list = [x for x in self.seeds if x.time>=start_date and x.time<=end_date]
        data_percentage=len(seed_list)/float(duration.days)*100.
        self.data_perentage=data_percentage
   
        for n in range(0,len(datelist)):
            if datelist[n] not in time_string_list:
                missing_date.append(datelist[n])
        self.missing_date = missing_date    
        
        # output the missing date to file
        if output_to_file:
            f=open('{0:s}_{1:s}_missing_data.txt'.format(self.network_code, self.code),'w')
            f.write("Station start date: {0:%Y/%m/%d %H:%M:%S.%f}\n".format(station.start_date.datetime))
            f.write('Station end date: {0:%Y/%m/%d %H:%M:%S.%f}\n'.format(station.end_date.datetime))
            f.write("Start date: {0:%Y/%m/%d %H:%M:%S.%f}\n".format(start_date.datetime))
            f.write('End date: {0:%Y/%m/%d %H:%M:%S.%f}\n'.format(end_date.datetime))
            for tmpdate in missing_date:
                f.write('{0:%Y/%m/%d}\n'.format(tmpdate.datetime))  
            f.close()           
            
        # check the size of existing seed files
        for n in range(0,len(size_string_list)):
            if size_string_list[n]<=3.0:
                small_date.append(time_string_list[n])
                small_size.append(size_string_list[n])
        self.small_date = small_date 
        # output seed file that is to small in size
        if output_to_file:
            f1=open('{0:s}_{1:s}_small_data.txt'.format(self.network_code, self.code),'w')
            for n in range(0,len(small_date)):
                f1.write('{0:%Y/%m/%d} {1:4.1f} MB\n'.format(small_date[n].datetime,small_size[n]))
            f1.close()
        # get the yearly completeness
        self._get_yearly_completeness()
        return self
    
    def _get_yearly_completeness(self,start_year=None,end_year=None):
        """check the data completeness for each year following the AGS convension
        e.g., 
        stn	2010 2011 2012 2013   2014 2015 2016 2017		
        BR2	N/A	 N/A  N/A  0.167	 N/A	  N/A  N/A N/A
        """
        # prepare date string list
        time_string_list = list()
        seed_list = list()
        missing_date = list()
        if start_year is None:       
            start_year = 2006

        if end_year is None:
            end_year = 2017
            
        for n in range(0,len(self.seeds)):
            tmp_time_string = self.seeds[n].time
            time_string_list.append(tmp_time_string.datetime)   
        # dictionary variable
        result = OrderedDict()
        result["Station"] = self.code
        for year in np.arange(start_year, end_year+1):
            start_date = UTCDateTime(year,1,1)
            end_date= UTCDateTime(year,12,31)
            duration=(end_date.datetime-start_date.datetime).days+1
            seed_list = [x for x in time_string_list if x>=start_date and x<=end_date]
            num_of_seed = len(seed_list)
            percentage=float(num_of_seed)/float(duration)*100.
            result["{0}".format(year)]="{0:3.1f}%".format(percentage)
        self.data_completeness=result
        return self

    def get_statistics(self):
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
        """
        Change Log:
        Mar. 4, 2018, Y.C., read individual station xml instead of network ones
        """
        from obspy import read_inventory
        inventory_path = parfile.stationxml_directory+self.network_code+'.'+self.code+'.xml'
        inv = read_inventory(inventory_path,format='STATIONXML')
        network = inv.networks[0].select(station=self.code)
        station = network[0]
        self.dataless = station
        self.latitude = station.latitude
        self.longitude = station.longitude
        self.elevation = station.elevation
        return self
    
    def select(self,time=None):
        """ Select seed that matches the given time, modified after 
        obspy select function
        ---------------------------------------------------------------------------
        Change Log:
        Feb. 28, 2018, Y.C., check if input date is UTCDateTime object
        ---------------------------------------------------------------------------  
        
        """                
        if not isinstance(time, UTCDateTime):
            warnings.warn('The input datestr is not a UTCDateTime object')
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
               "\tNetwork Code: {network_code}\n"
               "\tStation Code: {station_code}\n"
               "\tSeed Count: A total of {total} seeds\n"
               "\tStart-End: {start_date} - {end_date}\n"
               "\tLatitude: {lat:.2f}, Longitude: {lng:.2f}, "
               "Elevation: {elevation:.1f} m\n"
               "\tDataless: {dataless}\n")
        ret = ret.format(
            station_name=contents["stations"][0],
            network_code=self.network_code,
            station_code=self.code,
            total=self.total_number_of_seeds,
            start_date=str(self.start_date) if self.start_date else "",
            end_date=str(self.end_date) if self.end_date else "",
            lat=self.latitude if self.latitude else -999, 
            lng=self.longitude if self.longitude else -999,
            elevation=self.elevation if self.elevation else -999,
            dataless=self.dataless if self.dataless else '')
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
    def __init__(self,code='',path='',time='',size=''):
        self.code = code
        self.path = path
        self.time = time
        self.size = size
        
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
    
    def merge_seed(self,filein=None,target_directory=None,target_seed=None):
        """ Merge a daily seed/miniseed files with an target seed
        Change log:  Mar. 23, 2018, Y.C., add filein option 
        """
        import subprocess
        if not os.path.exists(self.path):  
            print "seed file doesn't exists"  
            raise IOError  
        if filein is None:
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
                            channel=None,if_merge=False):
        """ Convert the Fullseed file to miniseed
        could append argument to extract certain station or component
          -C arg retrieve the comments where 'arg' is either STN or CHN
        Change log: Mar. 23, 2018, Y.C., add merge option  
        """    
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
            if if_merge: # merge the file do not overwrite
                print new_filename
                self.merge_seed(filein=default_filename,
                                target_directory=output_directory,
                                target_seed=output_name)
                # delete miniseed
                os.unlink(default_filename)
            else: # overwrite the file
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
               "\tPath: {path}\n"
               "\tSize: {size:5.1f} MB\n")
        ret = ret.format(
            seed_code=self.code,
            date=str(self.time) if self.time else "",
            path=self.path,
            size=float(self.size)/1024./1024.)
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
