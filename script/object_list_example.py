#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 12:19:44 2017

@author: yunfeng
"""

class People(object):
    def __init__(self, people=()):
        self._members = list(people)

    # add other methods/special methods like __len__ etc.

    def add(self, person):
        self._members.append(person)
    def remove(self, person):
        self._members.remove(person)


people = People()
people.add('john')
people.add('jack')
