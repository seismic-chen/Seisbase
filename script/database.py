#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 09:24:07 2017

@author: yunfeng
"""
class Database(object):
    """ Database class contains Station class """

    def __init__(self, station):
        """ Create a new station class """
        self.station = station
        
class Station(object):
    """ Station class contains the Seed class """
    def __init__(self,seed=list(),code=''):
        self.seed = seed
        self.code = code
        
#    def get_statics(self):
#        """ Return data statistics for a given station """
#        self.
        
    def add(self,seedfile):
        """ add new Seed class """
        self.seed.append(seedfile)
        
class Seed(object):
    """ Seed class contains the FullSeed file information """
    def __init__(self,code='',path='',time=''):
        self.code = code
        self.path = path
        self.time = time
#    def distance_from_orogin(self):
#        """ computer my distance from the origin """
#        return ((self.x ** 2) + (self.y ** 2)) ** 0.5
#    def __add__(self, other):
#        new = copy.deepcopy(self)
#        new += other
#        return new
