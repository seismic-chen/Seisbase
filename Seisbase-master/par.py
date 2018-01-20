#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 19:18:07 2018
Par file that defines the path of the database structure
@author: yunfeng
"""
class Parfile(object):
    
    def __init__(self,database_directory='',inventory_path='',rdseed_path=''):
        self.database_directory=database_directory
        self.inventory_path=inventory_path
        self.rdseed_path=rdseed_path
        if not self.database_directory:
            self.database_directory='/Users/yunfeng/20_30/research/python_codes/'+\
                'seisbase_python/database'
        if not self.inventory_path:
            self.inventory_path='/Users/yunfeng/20_30/research/python_codes/'+\
                'seisbase_python/Seisbase-master/CRANE.xml'
        if not self.rdseed_path:
            self.rdseed_path='/usr/local/bin/rdseed'
#    def get_path(self):
#        if self.database_directory:
#            database_directory='/Users/yunfeng/20_30/research/python_codes/'+\
#                'seisbase_python/database'
#        if self.inventory_path:
#            inventory_path='/Users/yunfeng/20_30/research/python_codes/'+\
#                'seisbase_python/Seisbase-master/CRANE.xml'
#        return self
        
