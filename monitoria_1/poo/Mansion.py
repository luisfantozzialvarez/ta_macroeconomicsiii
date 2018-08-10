#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 19:09:22 2018

@author: luisfantozzialvarez
"""
from House import House

class Mansion(House):
    
    #Class contructor
    def __init__(self, sqr_metres, no_bathrooms, no_floors, no_rooms):
        #Calls constructor from parent. A mansion always has a pool!
        super().__init__(sqr_metres, no_bathrooms, no_floors, no_rooms, True)
        
    def create_description(self): 
        return super().create_description() + " This is a mansion!!!"
        