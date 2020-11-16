# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 15:38:05 2020

@author: Viktor Cheng
"""

import time
import pickle
import numpy as np
from numba import njit, jitclass, prange, boolean, int32, double
import math
import numpy as np
import scipy as sci
from scipy.optimize import minimize

from ruspy.estimation.estimation import estimate
from ruspy.simulation.simulation import simulate
from ruspy.model_code.fix_point_alg import calc_fixp
from ruspy.model_code.cost_functions import lin_cost
from ruspy.model_code.cost_functions import calc_obs_costs
from ruspy.estimation.estimation_transitions import create_transition_matrix




####Model#####

#to do:finish up model construction

class rust1987():
    
        def __init__(self,sol_gp=False,**kwargs): # called when created
        
        
        
        
        
        
        
        def solve(self, df=0.9999, s=90, cs=1e-3):
            """
            df stands for dicount factor
            s stands for states
            cs stands for cost of states
            """
            init_dict_nfxp = {
                "model_specifications": {
                "discount_factor": df,
                "number_states": s,
                "maint_cost_func": "linear",
                "cost_scale": cs,
                },
                "optimizer": {
                "approach": "NFXP",
                "algorithm": "scipy_L-BFGS-B",
                "gradient": "Yes",
                },
        
                }
            trans_result, fixed_result = estimate(init_dict_nfxp, model)
        
        def simulation(self):


            # Set simulating variables
            disc_fac = 0.9999
            num_buses = 100
            num_periods = 80000
            gridsize = 1000
            # We use the cost parameters and transition probabilities from the replication
            params = np.array([10.07780762, 2.29417622])
            trans_probs = np.array([0.39189182, 0.59529371, 0.01281447]) #maybe we can use the solve input here
            scale = 1e-3


            init_dict = {
            "simulation": {
            "discount_factor": disc_fac,
            "periods": num_periods,
            "seed": 123,
            "buses": num_buses,
                          },
             "plot": {"gridsize": gridsize},
            }

            num_states = 200 
            costs = calc_obs_costs(num_states, lin_cost, params, scale)

            trans_mat = create_transition_matrix(num_states, trans_probs)
            ev = calc_fixp(trans_mat, costs, disc_fac)[0]

            df = simulate(init_dict["simulation"], ev, costs, trans_mat)











#code waiting for integration:



class SimulatedMinimumDistance():
    ''' 
    This class performs simulated minimum distance (self) estimation.
    Input requirements are
    - model: Class with solution and simulation capabilities: model.solve() and model.simulate(). 
             Properties of model should be contained in model.par
    - mom_data: np.array (1d) of moments in the data to be used for estimation
    - mom_fun: function used to calculate moments in simulated data. Should return a 1d np.array
    
    '''    

    def __init__(self,model,mom_data,mom_fun,name='baseline',method='nelder-mead',est_par=[],lb=[],ub=[],options={'disp': False},print_iter=False,**kwargs): # called when created
        
        self.model = model
        self.mom_data = mom_data
        self.mom_fun = mom_fun
        self.name = name

        # estimation settings
        self.options = options
        self.print_iter = print_iter
        self.method = method

        self.lb = lb
        self.ub = ub

        self.est_par = est_par


    def obj_fun(self,theta,W,*args):
        
        if self.print_iter:
            for p in range(len(theta)):
                print(f' {self.est_par[p]}={theta[p]:2.3f}', end='')

        # 1. update parameters 
        for i in range(len(self.est_par)):
            setattr(self.model.par,self.est_par[i],theta[i]) # like par.key = val

        # 2. solve model with current parameters
        self.model.solve()

        # 3. simulate data from the model and calculate moments [have this as a complete function, used for standard errors]
        self.model.simulate()
        self.mom_sim = self.mom_fun(self.model.sim,*args)

        # 4. calculate objective function and return it
        diff = self.mom_data - self.mom_sim
        self.obj  = (np.transpose(diff) @ W) @ diff

        if self.print_iter:
            print(f' -> {self.obj:2.4f}')

        return self.obj 

    def estimate(self,theta0,W,*args):
        # TODO: consider multistart-loop with several algortihms - that could alternatively be hard-coded outside
        assert(len(W[0])==len(self.mom_data)) # check dimensions of W and mom_data

        # estimate
        self.est_out = minimize(self.obj_fun, theta0, (W, *args), method=self.method,options=self.options)

        # return output
        self.est = self.est_out.x
        self.W = W

    def std_error(self,theta,W,Omega,Nobs,Nsim,step=1.0e-4,*args):
        ''' Calculate standard errors and sensitivity measures '''

        num_par = len(theta)
        num_mom = len(W[0])

        # 1. numerical gradient. The objective function is (data - sim)'*W*(data - sim) so take the negative of mom_sim
        grad = np.empty((num_mom,num_par))
        for p in range(num_par):
            theta_now = theta[:] 

            step_now  = np.zeros(num_par)
            step_now[p] = np.fmax(step,step*theta_now[p])

            self.obj_fun(theta_now + step_now,W,*args)
            mom_forward = - self.mom_sim

            self.obj_fun(theta_now - step_now,W,*args)
            mom_backward = - self.mom_sim

            grad[:,p] = (mom_forward - mom_backward)/(2.0*step_now[p])

        # 1.1 reset the parameters in the model to theta
        for i in range(len(self.est_par)):
            setattr(self.model.par,self.est_par[i],theta[i]) 

        # 2. asymptotic standard errors [using Omega: V(mom_data_i). If bootstrapped, remember to multiply by Nobs]
        GW  = np.transpose(grad) @ W
        GWG = GW @ grad

        Avar = np.linalg.inv(GWG) @ ( GW @ Omega @ np.transpose(GW) ) @ np.linalg.inv(GWG)
        fac  = (1.0 + 1.0/Nsim)/Nobs # Nsim: number of simulated observations, Nobs: number of observations in data
        self.std = np.sqrt( fac*np.diag(Avar) )

        # 3. Sensitivity measures
        self.sens1 = - np.linalg.inv(GWG) @ GW  # Andrews I, Gentzkow M, Shapiro JM: "Measuring the Sensitivity of Parameter Estimates to Estimation Moments." Quarterly Journal of Economics. 2017;132 (4) :1553-1592
       

    def sensitivity(self,theta,W,fixed_par_str=None,step=1.0e-4,*args):
        ''' sensitivity measures '''

        num_par = len(theta)
        num_mom = len(W[0])

        # 1. numerical gradient. The objective function is (data - sim)'*W*(data - sim) so take the negative of mom_sim
        grad = np.empty((num_mom,num_par))
        for p in range(num_par):
            theta_now = theta[:] 

            step_now    = np.zeros(num_par)
            step_now[p] = np.fmax(step,step*theta_now[p])

            self.obj_fun(theta_now + step_now,W,*args)
            mom_forward = - self.mom_sim

            self.obj_fun(theta_now - step_now,W,*args)
            mom_backward = - self.mom_sim

            grad[:,p] = (mom_forward - mom_backward)/(2.0*step_now[p])

        # 1.1 reset the parameters in the model to theta
        #for i in range(len(self.est_par)):
        #    setattr(self.model.par,self.est_par[i],theta[i]) 

        # 2. Sensitivity measures
        GW  = np.transpose(grad) @ W
        GWG = GW @ grad
        Lambda = - np.linalg.inv(GWG) @ GW

        # 3. Sensitivity measures
        self.sens1 = Lambda  # Andrews I, Gentzkow M, Shapiro JM: "Measuring the Sensitivity of Parameter Estimates to Estimation Moments." Quarterly Journal of Economics. 2017;132 (4) :1553-1592

        # DO my suggestion
        if fixed_par_str:
            # mine: calculate the numerical gradient wrt parameters in fixed_par
            # update parameters
            for p in range(len(self.est_par)):
                setattr(self.model.par,self.est_par[p],theta[p])

            # change the estimation parameters to be the fixed ones
            est_par = self.est_par
            self.est_par = fixed_par_str

            # construct vector of fixed values
            gamma = np.empty(len(self.est_par))
            for p in range(len(self.est_par)):
                gamma[p] = getattr(self.model.par,self.est_par[p])

            # calculate gradient with respect to gamma
            num_gamma = len(gamma)
            grad_g = np.empty((num_mom,num_gamma))
            for p in range(num_gamma):
                gamma_now = gamma[:] 

                step_now    = np.zeros(num_gamma)
                step_now[p] = np.fmax(step,step*gamma_now[p])

                self.obj_fun(gamma_now + step_now,W,*args)
                mom_forward = - self.mom_sim

                self.obj_fun(gamma_now - step_now,W,*args)
                mom_backward = - self.mom_sim

                grad_g[:,p] = (mom_forward - mom_backward)/(2.0*step_now[p])

            # 1.1 reset the parameters in the model to theta
            for i in range(len(self.est_par)):
                setattr(self.model.par,self.est_par[i],gamma[i]) 

            self.est_par = est_par
            self.sens2 = Lambda @ grad_g

            ela = np.empty((len(theta),len(gamma)))
            for t in range(len(theta)):
                for g in range(len(gamma)):
                    ela[t,g] = self.sens2[t,g]*gamma[g]/theta[t]    
            
            self.sens2e = ela
