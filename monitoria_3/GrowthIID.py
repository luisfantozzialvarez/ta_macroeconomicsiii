#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:29:09 2018

@author: thiagoyuditachibana
"""
import numpy as np

class GrowthIID:
    
    
    #Class constructor: variables are:
    #eta - CRRA parameter
    #beta - discount factor
    #alpha - concavity parameter of production function
    #delta - depreciation rate
    #probz - vector with shock probabilities
    #valz - vector with shock realisation values
    #kgrid - grid with possible choices of capital
    
    def __init__(self, eta, beta, alpha, delta, probz, valz, kgrid):
        self.eta, self.beta, self.alpha, self.delta, self.probz, self.valz, self.kgrid =\
        eta, beta, alpha, delta, probz, valz, kgrid
        
        #Computes expected value of shock
        z_bar = np.dot(self.probz, self.valz)
        
        #Constructs dim Kprime x dim K x dim Z  K-array by replicating kgrid across columns and matrices
        #also constructs dim Kprime x dim K x dim Z Z-matrix by replicating zvalz across rows and matrices
        self.gridKprime, self.gridK, self.gridZ   = np.meshgrid(self.kgrid, self.kgrid, self.valz, indexing = 'ij')
        
        #Steady state of nonstochastic model with shock at average
        self.kss = ((z_bar*self.alpha)/((1/self.beta) + self.delta - 1))**(1/(1-self.alpha))
    
    #Function to reconstruct grids given a new vector of krigd
    def reconstr_grid(self, kgrid):
        self.kgrid = kgrid
        #Constructs dim Kprime x dim K x dim Z  K-array by replicating kgrid across columns and matrices
        #also constructs dim Kprime x dim K x dim Z Z-matrix by replicating zvalz across rows and matrices
        self.gridKprime, self.gridK, self.gridZ   = np.meshgrid(self.kgrid, self.kgrid, self.valz, indexing = 'ij')
    
    def utility(self, c):
        return c**(1- self.eta)/(1-self.eta)
    
    def prodf(self, k):
        return k**self.alpha
    
    
    #Bellman operator: given a dim K x dim Z matrix V(K,Z), computes TV
    def bellman(self, V):
        #Computes EV(K',Z')=
        EV = V @ self.probz
        
        #Reshapes EV to a cube
        EV_reshaped = np.repeat(EV, len(self.kgrid)*len(self.valz)).reshape((len(self.kgrid), len(self.kgrid), len(self.valz) ))
        
        #Computes consumption
        c = (1-self.delta)*self.gridK + self.gridZ*self.prodf(self.gridK) \
        - self.gridKprime
        
        #Sets consumption close to zero if negative
        c[c<=0] = 1e-10  
        
        #Computes utility
        uu = self.utility(c)  + self.beta*EV_reshaped
        
        return [np.max(uu, axis = 0),uu]
    
    #Runs Bellman and also returns policy function 
    #Value function iteration: Initial guess is a dimK x dimZ matrix V(K,z)
    def value_iteration(self, V0, tol = 1e-3):
        err = tol + 1
        while err > tol:
            Vnew, uu = self.bellman(V0)
            err = np.max(np.abs(Vnew - V0))
            V0 = Vnew
            print(err)
        
        policy_index = np.argmax(uu, axis = 0)
        policy_kprime = self.kgrid[policy_index]
        
        return [V0, policy_kprime, policy_index]
    
    def policy_iteration(self, V0, tol = 1e-3, Nh=10):
        #2D grids of (K,Z)
        gridK2D, gridZ2D = np.meshgrid(self.kgrid, self.valz, indexing = 'ij')
        err = tol + 1
        while err > tol:
            Vnew, uu = self.bellman(V0)
            policy_index = np.argmax(uu, axis = 0)
            
            for j in range(Nh):
                c =(1-self.delta)*gridK2D + gridZ2D*self.prodf(gridK2D) \
        - self.kgrid[policy_index]
                Vnew = self.utility(c) + self.beta*(Vnew@self.probz)[policy_index]
                
            err = np.max(np.abs(Vnew - V0))
            V0 = Vnew
            print(err)
        
        V0, uu =  self.bellman(V0)
        
        policy_index = np.argmax(uu, axis = 0)
        policy_kprime = self.kgrid[policy_index]
        
        return [V0, policy_kprime, policy_index]
        
            
            
            
        
    
        
    
    
        
        

