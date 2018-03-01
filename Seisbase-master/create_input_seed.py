"""
python code to create the input data
Nov 13, 2016, Yunfeng Chen, Global Seismology Group, University of Alberta
Nov. 14, 2017, Y.C., GSG, UofA.
Use database class and change the code to reflect more general use
Jan. 4, 2018, Y.C., merge dataless and miniseed files
Jan. 21, 2018, Y.C., modify the code to use new database class, which includes
network class
Feb. 24, 2018, Y.C., combine the data for any given number of stations
 
This code prepare miniseed data for ancc code, adopt resample routine from 
seismic ambient noise tomography package
"""
from obspy import read
import glob
import os
from database import Database,Station,Seed
import pickle 
from obspy.io.xseed import Parser
import numpy as np
from par import Parfile
parfile=Parfile()
#import resample
def resample(trace, dt_resample):
    """
    Subroutine to resample trace from seismic ambient noise package
    @type trace: L{obspy.core.trace.Trace}
    @type dt_resample: float
    @rtype: L{obspy.core.trace.Trace}
    """
    dt = 1.0 / trace.stats.sampling_rate
    factor = dt_resample / dt
    if int(factor) == factor:
        # simple decimation (no filt because it shifts the data)
        trace.decimate(int(factor), no_filter=True)
    else:
        # linear interpolation
        tp = np.arange(0, trace.stats.npts) * trace.stats.delta
        zp = trace.data
        ninterp = int(max(tp) / dt_resample) + 1
        tinterp = np.arange(0, ninterp) * dt_resample

        trace.data = np.interp(tinterp, tp, zp)
        trace.stats.npts = ninterp
        trace.stats.delta = dt_resample
        trace.stats.sampling_rate = 1.0 / dt_resample
        #trace.stats.endtime = trace.stats.endtime + max(tinterp)-max(tp)

# load the database object
file_db = open('noise_database.obj', 'r') 
db = pickle.load(file_db)

# define a station/network list that need to be cross-correlated with
stalist={'ANMO','DWPF','COR'}
netlist={'IU','IU','IU'}


# define two stations that need to be cross-correlated with
sta1='ANMO'
net1='IU'
sta2='DWPF'
net2='IU'
db_sta1 = db.select(network=net1,station=sta1)
db_sta2 = db.select(network=net2,station=sta2)


miniseed_file='{0:s}_{1:s}.seed'.format(sta1,sta2)
# convert the station xml to dateless seed
stationxml_directory=parfile.stationxml_directory

xml1 = stationxml_directory+net1+'.'+sta1+'.xml'
xml2 = stationxml_directory+net2+'.'+sta2+'.xml'
convertor_path=stationxml_directory+'stationxml-converter.jar'
dataless_file=sta1+'_'+sta2+'.dataless'
cmd='java -jar {0:s} --seed {1:s} {2:s} --output {3:s}'.format(convertor_path,xml1,xml2,dataless_file)
#    sp1 = Parser('/Users/yunfeng/20_30/research/python_codes/seismic-noise-tomography-master/'+
#                'StationXML/'+st1[0].stats['station']+'.dataless')
#    sp2 = Parser('/Users/yunfeng/20_30/research/python_codes/seismic-noise-tomography-master/'+
#                'StationXML/'+st2[0].stats['station']+'.dataless')
os.system(cmd)
# filter options
freqmin=1./200.
freqmax=1./3.
corners=2
zerophase=True
period_resample=1.0
#for BRU_data in glob.glob(database_directory + 'BRU/2006/*seed'):
for seed1 in db_sta1.networks[0].stations[0].seeds[0:365]: 
    # check if the corresponding data from station2 exists
    seed2_list = db_sta2.networks[0].stations[0].select(time=seed1.time).seeds
    if seed2_list:  # check if list is empty
        seed2=seed2_list[0]
        if os.path.isfile(seed2.path):
            print str(seed2.time) + "Seed file from station 2 exists"
            try:
                year='{:04d}'.format(seed1.time.year)
                month='{:d}'.format(seed1.time.month)
                output_directory = '/Users/yunfeng/20_30/research/ambient_noise/'+\
                'ancc-1.0-0.src/SEED/'+year+'/'+month
                # create directory if not exits
                if not os.path.exists(output_directory):
                    os.makedirs(output_directory)
                st1 = []
                st1 = seed1.read_seed()
            except:
                pass
            # check if successfully read
            if not st1:
                continue
            # change the network code
            st2 = []
            try:
                st2 = seed2.read_seed()
            except:
                pass
            # check if successfully read
            if not st2:
                continue	
            # resample the data to 1Hz
            # follow the pscrosscorr precedure, filter first then resample
            st1.filter(type="bandpass",
                       freqmin=freqmin,
                       freqmax=freqmax,
                       corners=corners,
                       zerophase=zerophase)
            st2.filter(type="bandpass",
                       freqmin=freqmin,
                       freqmax=freqmax,
                       corners=corners,
                       zerophase=zerophase)
            for n1 in range(0,len(st1)):
                resample(st1[n1], dt_resample=period_resample)
                # change the network name
                if st1[n1].stats.network == 'XX':
                    st1[n1].stats.network = 'CR'
            for n2 in range(0,len(st2)):
                resample(st2[n2], dt_resample=period_resample)
                if st2[n2].stats.network == 'XX':
                    st2[n2].stats.network = 'CR'
#    sampling_rate = int(st1[0].stats.sampling_rate)
#    st1.decimate(sampling_rate,strict_length=False)
#    sampling_rate = int(st2[0].stats.sampling_rate)
#    st2.decimate(sampling_rate,strict_length=False)
    # combine seed files from these two stations
        for n in range(0,len(st2)):
            st1.append(st2[n])
        st1.write(miniseed_file,format='MSEED')
        # convert the miniseed file to full seed
        cmd1='/usr/local/bin/rdseed -d -o 5 -f {0:s} -g {1:s} -l >> output.txt'.format(miniseed_file,dataless_file)
        os.system(cmd1)
        # rename the file
        year = st1[0].stats.starttime.year
        month = st1[0].stats.starttime.month
        day = st1[0].stats.starttime.day
        jday = st1[0].stats.starttime.julday
        seed_file = 'D.{0:4d}.{1:03d}.{2:d}.{3:d}.000000.seed'.format(year,jday,month,day)
        # check if seed file exists
        if os.path.isfile('/Users/yunfeng/20_30/research/python_codes/seisbase_python/Seisbase-master/seed.rdseed'):
            os.rename('seed.rdseed',output_directory+'/'+seed_file)