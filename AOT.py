# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 11:52:52 2021

@author: Amy
"""

"""define parameters"""
S_0 = 100
K = 100
r = 0.015
T = 0.5
D = 120

"""import module"""
import math
import numpy as np
    
class uoc_BlackScholes:
    """
    Class containing data and methods to describe 
    the Black−Scholes PDE for up-and-out call option
    """
    def __init__(self, rate, strike, barrier, time, alpha, S_low, S_up):
        """ Constructor
        : param rate: risk-free rate
        : param strike: option’s strike price
        : param barrier: option’s up-and-out barrier
        : param time: option’s time to expiration
        : param alpha: use to compute standard deviation of option
        : param S_low: lower bound for the option's spot price
        : param S_up: upper bound for the option’s spot price
        """
        self.rate = rate
        self.strike = strike
        self.barrier = barrier
        self.time = time
        self.alpha = alpha
        self.S_low = S_low
        self.S_up = S_up
        
    def sigma(self, t, S):
        """ Compute standard deviation
        : param t: current time 
        : param S: current spot price
        """
        return 0.13*math.exp(-t)*((100/S)**self.alpha)
    
    def coeff_a(self, t, S):
        """ Compute coefficient a
        : param t: current time 
        : param S: current spot price
        """
        return -(((self.sigma(t, S))**2)*(S**2))/2
    
    def coeff_b(self, t, S):
        """ Compute coefficient b
        : param t: current time 
        : param S: current spot price
        """
        return -self.rate*S
    
    def coeff_c(self, t, S) :
        """ Compute coefficient c
        : param t: current time 
        : param S: current spot price
        """
        return self.rate
    
    def coeff_d(self, t, S) :
        """ Compute coefficient d
        : param t: current time 
        : param S: current spot price
        """
        return 0
    
    def bound_payoff(self, S) :
        """ Compute boundary condition for payoff
        : param S: current spot price
        """
        if S < self.barrier:
            return max(S - self.strike, 0)
        else:
            return 0
    
    def bound_Slow(self, t) :
        """ Compute payoff for lowest spot price
        """
        return 0
    
    def bound_Sup(self, t) :
        """ Compute payoff for highest spot price
        """
        return 0

"""define PDE that is needed to be solved in (a)"""
PDE_a = uoc_BlackScholes(
    rate = r, 
    strike = K, 
    barrier = D, 
    time = T, 
    alpha = 0.2, 
    S_low = 0, 
    S_up = D
    )

class PDESolver :
    """Abstract class to solve Black−Scholes PDEs."""
    def __init__(self, pde, imax, jmax):
        """Constructor
        : param pde: The PDE to solve
        : param imax: last value of the first variable’s discretisation
        : param jmax: last value of the second variable’s discretisation
        """
        self.pde = pde
        self.imax = imax
        self.jmax = jmax
        self.dt = self.pde.time / imax
        self.dS = (self.pde.S_up-self.pde.S_low) / jmax
        """grid is used to save solutions later"""
        self.grid = np.empty((self.imax+1, self.jmax+1),dtype=float)
        
    def t(self, i):
        """Return the discretised value of t at index i"""
        return self.dt * i
    
    def S(self, j):
        """Return the discretised value of S at index j"""
        return self.dS * j + self.pde.S_low
    
    def new_a(self, i, j):
        """
        Helper umbrella function to get coefficient a at discretised locations
        """
        return self.pde.coeff_a(self.t(i),self.S(j))
    
    def new_b(self, i, j):
        """
        Helper umbrella function to get coefficient b at discretised locations
        """
        return self.pde.coeff_b(self.t(i),self.S(j))
    
    def new_c(self, i, j):
        """
        Helper umbrella function to get coefficient c at discretised locations
        """
        return self.pde.coeff_c(self.t(i),self.S(j))
    
    def new_d(self, i, j):
        """
        Helper umbrella function to get coefficient d at discretised locations
        """
        return self.pde.coeff_d(self.t(i),self.S(j))
    
    def new_payoff(self, j):
        """
        Helper umbrella function to get payoff boundary condition
        for time at discretised j index
        : param j: discretised S index
        """
        return self.pde.bound_payoff(self.S(j))
    
    def new_Slow(self, i):
        """
        Helper umbrella function to get lower boundary condition on spot price
        for time at discretised i index
        : param j: discretised t index
        """
        return self.pde.bound_Slow(self.t(i))
    
    def new_Sup(self, i):
        """
        Helper umbrella function to get upper boundary condition on spot price
        for time at discretised i index
        : param j: discretised t index
        """
        return self.pde.bound_Sup(self.t(i))
    
    def interpolate(self, t, S):
        """
        Get interpolated solution value at given time and space
        : param t: point in time
        : param S: point in space
        : return: interpolated solution value
        """
        i = int(t / self.dt)
        j = int((S - self.pde.S_low) / self.dS)
        x1 = t / self.dt - i
        x0 = 1 -x1
        y1 = (S - (self.pde.S_low + self.dS * j)) /self.dS
        y0 = 1 - y1
        return (x1 * y1 * self.grid[i + 1,j + 1]
                + x1 * y0 * self.grid[i + 1,j]
                + x0 * y1 * self.grid[i,j + 1]
                + x0 * y0 * self.grid[i,j])
    
class Explicit_Scheme(PDESolver):
    """ Black−Scholes PDE solver using the explicit scheme
    """
    def __init__(self, pde, imax, jmax):
        super().__init__(pde, imax, jmax)
        
    def coeff_A(self, i, j):
        """
        Coefficient A {i, j} for the explicit scheme
        : param i : index of S distretisation
        : param j : index of t discretisation
        """
        return self.dt/self.dS * (self.new_b(i,j)/2 - self.new_a(i,j)/self.dS)
    
    def coeff_B(self, i, j):
        """
        Coefficient B {i, j} for the explicit scheme
        : param i : index of S distretisation
        : param j : index of t discretisation
        """
        return (1 - self.dt*self.new_c(i, j) 
                + (2 * self.dt * self.new_a(i,j)) / (self.dS ** 2))

    def coeff_C(self, i, j):
        """
        Coefficient C {i, j} for the explicit scheme
        : param i : index of S distretisation
        : param j : index of t discretisation
        """
        return -self.dt/self.dS * (self.new_b(i,j)/2 + self.new_a(i,j)/self.dS)

    def coeff_D(self, i, j):
        """
        Coefficient D {i, j} for the explicit scheme
        : param i : index of S distretisation
        : param j : index of t discretisation
        """
        return -self.dt * self.new_d(i, j)

    def solve_grid(self):
        """
        Solves the PDE and saves the values in the matrix ‘grid‘
        """
        def update(i, j):
            """ Compute udpate for iterations"""
            return (self.coeff_A(i, j) * self.grid[i,j - 1]
                    + self.coeff_B(i, j) * self.grid[i,j]
                    + self.coeff_C(i, j) * self.grid[i,j + 1]
                    + self.coeff_D(i, j))
        
        self.grid[self.imax] = [self.new_payoff(j) for j in range(self.jmax+1)]
        
        for i in range(self.imax, -1, -1):
            self.grid[i-1,0] = self.new_Slow(i-1)
            self.grid[i-1,self.jmax] = self.new_Sup(i-1)
            self.grid[i-1,1:-1] = [update(i,j) for j in range(1,self.jmax)]


"""Test imax,jmax give same result"""
exp_scheme1 = Explicit_Scheme(PDE_a, 200, 200)
exp_scheme1.solve_grid()
soln_exp1 = exp_scheme1.interpolate(0, S_0)
print("Solution from the explicit scheme : {0}". format(soln_exp1))

exp_scheme2 = Explicit_Scheme(PDE_a, 90, 90)
exp_scheme2.solve_grid()
soln_exp2 = exp_scheme2.interpolate(0, S_0)
print("Solution from the explicit scheme : {0}". format(soln_exp2))

exp_scheme3 = Explicit_Scheme(PDE_a, 100, 100)
exp_scheme3.solve_grid()
soln_exp3 = exp_scheme3.interpolate(0, S_0)
print("Solution from the explicit scheme : {0}". format(soln_exp3))

exp_scheme4 = Explicit_Scheme(PDE_a, 150, 150)
exp_scheme4.solve_grid()
soln_exp4 = exp_scheme4.interpolate(0, S_0)
print("Solution from the explicit scheme : {0}". format(soln_exp4))




        
        





