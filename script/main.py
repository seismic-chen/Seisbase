#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 23:54:50 2017

@author: yunfeng
"""
#from point import Point
#p = Point(3,4)
#q = Point(5,6)
#print p.x
#print q.x
#print q.distance_from_orogin()
import os
import fnmatch
import re
from obspy import UTCDateTime
from database import Database,Station,Seed
st = Station()
database_directory='/Users/yunfeng/20_30/study/GASS_project/seisbase_python/database'
pattern='FullSeed_*.seed'
for root, dirs, files in os.walk(database_directory):
    for filename in fnmatch.filter(files,pattern):
        path=os.path.join(root,filename)
        datestr=re.split('_|\.',filename)[2]
        timestr=re.split('_|\.',filename)[3]
        se=Seed(filename,path,UTCDateTime(datestr+timestr))
        st.seed.append(se)

#se1 = Seed()
#se1.code = 'FullSeed_JOF_20070109_000000.seed'
#se2 = Seed()
#se2.code = 'FullSeed_JOF_20070110_000000.seed'

#st.seed.append(se1)
#st.seed.append(se2)
db = Database(st)