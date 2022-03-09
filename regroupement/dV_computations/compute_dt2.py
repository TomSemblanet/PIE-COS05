#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 7 09:12:33 2022

@author: g.pierre
"""

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd
from utils import debris_data_loader as DDL
from utils import constants
from matplotlib.colors import LinearSegmentedColormap

def propagate_RAAN(RAAN_0,a,i, t):
	return (-3/2*np.sqrt(constants.mu_EARTH/a**3)*((constants.R_EARTH/a)**2)*constants.J2*np.cos(i)*t + RAAN_0) % (2*np.pi)

# def compute_global_dt(debris):

# 	'''
# 	Arguments:
# 		debris (array): array containing orbital parameters of the debris sequence 

# 	Returns: 
# 		dt (list of float): the dt associated with each change of RAAN
# 	'''

# 	dt=[]
# 	omega=np.zeros(debris.index)
# 	omega[:,0]=debris.values[:][3]
# 	current_omega_i=omega[0,0]
# 	current_omega_f=omega[1,0]

# 	for i in range(debris.index-1):
# 		a_i,i_i=debris.values[i][[0, 2]]
# 		a_f,i_f=debris.values[i+1][[0, 2]]

# 		for j in range(i):
			
# 			current_omega_i=propagate_RAAN(omega[i,i-1],a_i,i_i,dt[i-1])
# 			current_omega_f=propagate_RAAN(omega[i+1,i-1],a_i,i_i,dt[i-1])
# 			omega[i,i]=current_omega_i
# 			omega[i+1,i]=current_omega_f

# 		dt.append(compute_dt(a_i,a_f,i_i,i_f,current_omega_i,current_omega_f))

# 	return dt

def compute_and_propagate(debris_op):

	N=np.size(debris_op,0)

	omega=np.zeros((N,N))
	omega[:,0]=debris_op[:][3]


	for i in range (N-1):

		dt.append(compute_dt(debris_op[i][0], debris_op[j][2],current_omega_i,current_omega_f))

		for j in range (i+1,N-1):	

			dt.append(compute_dt(debris_op[i][0], debris_op[i+1][0],debris_op[i][2], debris_op[i+1][2],omega[i][i],omega[i+1][i]))
			omega[j, i+1]= propagate_RAAN(omega[j,i], debris_op[j][0], debris_op[j][2], dt[i])


return dt

