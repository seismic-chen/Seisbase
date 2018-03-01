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
import datetime
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
stalist=['ANMO','DWPF','COR']
netlist=['IU','IU','IU']


miniseed_file='all_stations.seed'
# convert the station xml to dateless seed
stationxml_directory=parfile.stationxml_directory

xmllist=[stationxml_directory+netlist[n]+'.'+stalist[n]+'.xml' for n in range(0,len(netlist))]

convertor_path=stationxml_directory+'stationxml-converter.jar'
dataless_file='all_stations.dataless'
cmd='java -jar {0:s} --seed {1:s} --output {2:s}'.format(convertor_path,' '.join(xmllist),dataless_file)

os.system(cmd)
# filter options
freqmin=1./200.
freqmax=1./3.
corners=2
zerophase=True
period_resample=1.0
#for BRU_data in glob.glob(database_directory + 'BRU/2006/*seed'):
basetime=datetime.datetime(2006, 3, 10, 0, 0, 0)
numdays=1
timelist=[basetime + datetime.timedelta(days=x) for x in range(0, numdays)]
for t in timelist:
    stlist=read()
    stlist.clear() 
    for n in range(0,len(stalist)):
        db_sta = db.select(network=netlist[n],station=stalist[n])
        seedlist = db_sta.networks[0].stations[0].select(time=t).seeds
        if seedlist:  # check if list is empty
            seed=seedlist[0]
            if os.path.isfile(seed.path):
                print str(seed.time) + "Seed file from station {0:s} exists".format(stalist[n])
                try:
                    year='{:04d}'.format(seed.time.year)
                    month='{:d}'.format(seed.time.month)
                    output_directory = '/Users/yunfeng/20_30/research/ambient_noise/'+\
                    'ancc-1.0-0.src/SEED/'+year+'/'+month
                    # create directory if not exits
                    if not os.path.exists(output_directory):
                        os.makedirs(output_directory)
                    st = []
                    st = seed.read_seed()
                except:
                    pass
                # check if successfully read
                if not st:
                    continue               
                # follow the pscrosscorr precedure, filter first then resample
                st.filter(type="bandpass",
                           freqmin=freqmin,
                           freqmax=freqmax,
                           corners=corners,
                           zerophase=zerophase)
                for n1 in range(0,len(st)):
                    resample(st[n1], dt_resample=period_resample)
                    # change the network name
                    if st[n1].stats.network == 'XX':
                        st[n1].stats.network = 'CR'
         # append the trace from current station to stream list
        for n1 in range(0,len(st)):
             stlist.append(st[n1])               
                

    stlist.write(miniseed_file,format='MSEED')
    # convert the miniseed file to full seed
    cmd1='/usr/local/bin/rdseed -d -o 5 -f {0:s} -g {1:s} -l >> output.txt'.format(miniseed_file,dataless_file)
    os.system(cmd1)
    # rename the file
    year = stlist[0].stats.starttime.year
    month = stlist[0].stats.starttime.month
    day = stlist[0].stats.starttime.day
    jday = stlist[0].stats.starttime.julday
    seed_file = 'D.{0:4d}.{1:03d}.{2:d}.{3:d}.000000.seed'.format(year,jday,month,day)
    # check if seed file exists
    if os.path.isfile('/Users/yunfeng/20_30/research/python_codes/seisbase_python/Seisbase-master/seed.rdseed'):
        os.rename('seed.rdseed',output_directory+'/'+seed_file)