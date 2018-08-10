#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 18:13:01 2018

@author: luisfantozzialvarez
"""
from House import House
from Mansion import Mansion
import numpy as np

luis_house = House(100, 4, 3, 4, True)
thiago_house = House(200, 10, 2, 3, False)

print(luis_house.create_description())
print(thiago_house.create_description())

bill_house = Mansion(1000, 3, 4, 5)
print(bill_house.size_cm2())
print(bill_house.create_description())