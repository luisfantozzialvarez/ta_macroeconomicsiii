'''
Created on 19 de ago de 2017

@author: luisfantozzialvarez
'''
import numpy as np
from scipy.optimize import fsolve

class Growth(object):

    #class constructor: initialises economy given parameters
    def __init__(self, beta, delta, gamma, eta, alpha, theta, tau_c, tau_h, tau_k, k0):
        self.beta, self.delta, self.gamma, self.eta, self.alpha, self.theta, self.tau_c, self.tau_h, self.tau_k, self.k0 = \
        beta, delta, gamma, eta, alpha, theta, tau_c, tau_h, tau_k, k0

    #the function below returns, in this order, the following variables in steady state:
    #h*, k*, c*, y*, i*, g* where y*, c*, k*, i*, g* are normalized 
    #per units of effective labour (At.Nt)
    def get_steady_state(self):
        #computes h*/k*
        h_div_k =  (((1+self.gamma)/(self.beta) + self.delta - 1)/((1- self.tau_k)*self.alpha))**(1/(1-self.alpha))
        #computes k*
        h_ss = 1 + (self.theta/((1- self.tau_h)*(1-self.alpha)))*((1-self.delta)*(h_div_k)**(self.alpha -1 ) - \
        (1+ self.eta)*(1+self.gamma)*(h_div_k)**(self.alpha-1) + (1 - self.tau_k)*self.alpha + (1-self.alpha)*(1-self.tau_h)) 
        h_ss = h_ss**(-1)
        #computes k_ss
        k_ss = h_ss/h_div_k
        #computes c_ss
        c_ss = ((1- h_ss)*(1-self.tau_h)*(1-self.alpha)*(h_div_k)**(-self.alpha))/(self.theta*(1+self.tau_c))
        #computes y_ss = k_ss^(alpha0*h_ss*(1-alpha)
        y_ss = k_ss**(self.alpha)*h_ss**(1-self.alpha)
        #computes i_ss
        i_ss = k_ss*((1+self.eta)*(1+self.gamma) - (1-self.delta))
        g_ss = self.tau_c*c_ss + self.tau_h*(1-self.alpha)*y_ss + self.tau_k*self.alpha*y_ss 
        output = np.array([h_ss, k_ss, c_ss, y_ss, i_ss, g_ss])
        return output
    
    #for given values of k_t, k_t+1, h_t and h_t+1, this returns the movement equations that characterise the economy 
    def eq_transition(self, k_t, k_t_1, h_t, h_t_1):
        output = np.empty(2)
        output[0] = ((1-h_t_1)/(1-h_t))*((k_t_1/h_t_1)**(self.alpha))*((k_t/h_t)**(-self.alpha)) \
        - (self.beta/(1+self.gamma))*(1 - self.delta + (1-self.tau_k)*self.alpha*(h_t_1/k_t_1)**(1-self.alpha))
        output[1] = ((1-h_t)*(1-self.tau_h)*(1-self.alpha)/self.theta)*(k_t/h_t)**self.alpha + (1+ self.eta)*(1+ self.gamma)*k_t_1 \
        -(1- self.delta + (1-self.tau_k)*self.alpha*(h_t/k_t)**(1-self.alpha))*k_t - (1-self.tau_h)*(1-self.alpha)*k_t**(self.alpha)*h_t**(1-self.alpha)
        return(output)
    
    #give me a vector x of dimension 2*T+2 and I will return
    #the movement equations of the economy from 0 to T, where k0 was given in the constructor,
    #and k_{t+1} is in its steady state value.
    #First T entries of x should be the path of capital. Remaining entries should be the path of h.
    def gen_transition(self, x):
        output = np.empty(len(x))
        T=int((len(x)-2)/2)
        for t in range(T+1):
            if t == 0:
                aux = self.eq_transition(self.k0, x[t], x[t+T], x[t+1+T])
                output[t] = aux[0]
                output[t+T+1] = aux[1]
            else:
                if t == T:
                    k_ss = self.get_steady_state()[1]
                    aux = self.eq_transition(x[t-1], k_ss, x[t+T], x[t+1+T])
                    output[t] = aux[0]
                    output[t+T+1] = aux[1]
                else:
                    aux = self.eq_transition(x[t-1], x[t], x[t+T], x[t+1+T])
                    output[t] = aux[0]
                    output[t+T+1] = aux[1]
        return(output)       
    
    #given paths for capital and hours worked, returns path of consumption per units of effective labour
    def gen_consumption_path(self, k,h):
        return ((1-h)*(1-self.tau_h)*(1-self.alpha)*(k/h)**self.alpha)/(self.theta*(1+self.tau_c))
    
    #given paths for capital and hours worked, returns path of output per units of effective labour
    def gen_output_path(self, k, h):
        return k**(self.alpha)*h**(1-self.alpha)
    
    #given path for capital, returns path of investment per units of effective labour. Last period investment 
    #is set to steady state value
    def gen_investment_path(self, k):
        i = k[1:len(k)]*(1+self.eta)*(1+self.gamma) - (1-self.delta)*k[0:len(k)-1] 
        i_ss = self.get_steady_state()[4]
        i = np.concatenate((i,[i_ss]))
        return i
    
    #given path for capital and hours worked, returns path of government spending per units of effective labour. 
    def gen_government_path(self, k_path, h_path):
        return self.tau_c*self.gen_consumption_path(k_path, h_path) + self.tau_h*(1-self.alpha)*self.gen_output_path(k_path, h_path) \
        + self.tau_k*self.alpha*self.gen_output_path(k_path, h_path)
    
    #solves model for T periods, where T+1 is steady state and t=0 is staring period
    #initial guess for k is a linear interpolation from k_0 to k_{steady\_state}
    #initial guess for h is constant at steady state
    #returns paths of the following variables, in the following order: h, k, c, y, i
    def solve_model(self, T):
        ss = self.get_steady_state()
        k_ss = ss[1]
        h_ss = ss[0]
        guess_k = np.linspace(self.k0, k_ss, T)
        guess_h = np.repeat(h_ss, T+2)
        guess = np.concatenate((guess_k, guess_h))
        
        solucao = fsolve(self.gen_transition, guess)
        
        k_path = solucao[0:T]
        k_path = np.concatenate(([self.k0], k_path, [k_ss]))
        
        h_path = solucao[T:]
        
        c_path = self.gen_consumption_path(k_path, h_path)
        y_path = self.gen_output_path(k_path, h_path)
        i_path = self.gen_investment_path(k_path)
        g_path = self.gen_government_path(k_path, h_path)
        return [h_path,k_path,c_path, y_path, i_path, g_path]
    
    
    #presents the solution of the model per units of N_t (labour), and not the default N_t.A_t (effective labour)
    #takes as arguments the number of periods where the approximation is made
    #and initial guess for A_0
    def solve_model_per_capita(self,T, A0):
        solution = self.solve_model(T)
        t =  np.linspace(0, T+1 , T+2)
        A = A0*(1+self.gamma)**t
        
        for i in  (range(len(solution)-1)):
            solution[i+1]= solution[i+1]*A
        return solution