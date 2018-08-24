#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 03:01:13 2018

"""
import numpy as np

#Meshgrid

k      = np.linspace(1,10,10)
kprime = np.linspace(1,10,10)
z      = np.array([1,2,3,4])

#(1) Value function format
kk, zz = np.meshgrid(k,z, indexing = 'ij')

#(2) Array for EV (to reduce computing time)

kkkprime, kkk, zzz = np.meshgrid(k,k,z, indexing = 'ij')