#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 19:29:46 2018

@author: luisfantozzialvarez
"""
import numpy as np
from scipy.optimize import fsolve

class DiscreteGrowthCapital:
    #Class consctructor: takes as arguments
    #A - productivity parameter
    #alpha - production function concavity
    #delta - depreciation rate
    #beta - discount factor
    #eta - risk-aversion
    #k0 - initial capital stock
    def __init__(self, A, alpha, delta, beta, eta, k0):
        self.A = A
        self.alpha = alpha
        self.delta = delta
        self.beta = beta
        self.eta = eta
        self.k0 = k0
        #Computes steady state capital stock
        self.kss = ((self.alpha*self.A)/((1/self.beta)- (1-self.delta)))**(1/(1-self.alpha))
    
    #Production function
    def f(self, k):
        return k**self.alpha
    
    #Marginal product
    def fprime(self, k):
        return self.alpha*k**(self.alpha-1)
    
    def compute_consumption(self, k, kprime):
        return self.f(k) + (1-self.delta)*k - kprime
    
    def uprime(self, c):
        return c**(-self.eta)
    
    #Computes the difference equations that determine equilibrium. Takes as arguments:
    #k.path - capital path from time 1 to time T (T+1 assumed to be steady state, k0 given)
    def system_equations(self, k_path):
        k_extended = np.concatenate(([self.k0],k_path, [self.kss]))
        k_current = k_path
        k_before = k_extended[:-2]
        k_after = k_extended[2:]
        equations = self.uprime(self.compute_consumption(k_before, k_current)) -\
        self.beta*(self.fprime(k_current) + (1-self.delta))*self.uprime(self.compute_consumption(k_current, k_after))
        return equations
    
    #Solves the model using the secant fsolve
    def solve_model(self, k_path):
        k_solution = fsolve(self.system_equations, k_path)
        return k_solution
        

    
