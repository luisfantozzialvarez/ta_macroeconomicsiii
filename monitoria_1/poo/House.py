#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 17:57:32 2018

@author: luisfantozzialvarez
"""

#this class creates a house
class House:  
            
    #Class variable -- this does NOT vary with instance
    what_is = "This is a house" 
    
    #Class contructor: this initialises the variables required by the object
    def __init__(self, sqr_metres, no_bathrooms, no_floors, no_rooms, pool):
        #Creates variables specific to this instance (object). This will vary
        #across objects
        self.sqr_metres = sqr_metres
        self.no_bathrooms = no_bathrooms
        self.no_floors = no_floors
        self.no_rooms = no_rooms
        self.pool = pool
        
    #Method that generates a description
    def create_description(self):
        result = "This house has " + str(self.sqr_metres) + " square metres, It has " + \
        str(self.no_bathrooms) + " bathrooms."
        
        if self.pool:
            result = result + "It has a swimming pool."
        else:
            result = result + "It does not have a swimming pool."
            
        return result
    
    
    def size_cm2(self):
        return self.sqr_metres*100**2

        

    
    


        
    