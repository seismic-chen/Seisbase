#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 22:30:52 2018
Class of noise correlation database
@author: yunfeng
"""
class Noisedata(object):
    """ Noise database class """
    
    def __init__(self, networks=list(),code=''):
        """ Create a new networks class """
        self.networks = networks
    def add(self,networkfile):
        """ add new network class """
        self.networks.append(networkfile)
