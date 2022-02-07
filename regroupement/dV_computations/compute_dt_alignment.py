#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 18:09:33 2022

@author: g.pierre
"""

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd
from utils import debris_data_loader as DDL
from utils import constants
from matplotlib.colors import LinearSegmentedColormap

def compute_dt(a1,a2,i1,i2,RAAN1, RAAN2, print_result=False):
    """ Function computing the delta_t [s] required to modify the RAAN of the orbit from RAAN1 to RAAN2 

            inputs : 
            ------
                    - a1 : float
                            Initial SMA [km]
                    - a2 : float
                            Final SMA [km]
                    - i1 : float
                            Initial inclination [rad]
                    - i2 : float
                            Final inclination [rad]
                    - RAAN1 : float
                            Initial Right ascension of the ascending node [rad]
                    - RAAN2 : float
                            Final Right ascension of the ascending node [rad]
            output : 
            ------
                    - dt : float
                            Required delta_t [s]
    """
    if ((RAAN1>RAAN2) and (constants.J2*(constants.R_EARTH**2)*((np.sqrt(constants.mu_EARTH/(a1**3))*np.cos(i1)/(a1**2))-(np.sqrt(constants.mu_EARTH/(a2**3))*np.cos(i2)/(a2**2))))>0 ):
        dt=-(2/3)*(RAAN1-RAAN2-2*np.pi)/(constants.J2*(constants.R_EARTH**2)*((np.sqrt(constants.mu_EARTH/(a1**3))*np.cos(i1)/(a1**2))-(np.sqrt(constants.mu_EARTH/(a2**3))*np.cos(i2)/(a2**2))))
        if print_result == True:
            print('DeltaT for RAAN1 = ', RAAN1, ' and RAAN2 = ', RAAN2, ' is ', dt/86400, ' days.')
        return dt
    elif ((RAAN2>RAAN1) and (constants.J2*(constants.R_EARTH**2)*((np.sqrt(constants.mu_EARTH/(a1**3))*np.cos(i1)/(a1**2))-(np.sqrt(constants.mu_EARTH/(a2**3))*np.cos(i2)/(a2**2))))<0): 
        dt= -(2/3)*(RAAN1+2*np.pi-RAAN2)/(constants.J2*(constants.R_EARTH**2)*((np.sqrt(constants.mu_EARTH/(a1**3))*np.cos(i1)/(a1**2))-(np.sqrt(constants.mu_EARTH/(a2**3))*np.cos(i2)/(a2**2))))
        if print_result == True:
            print('DeltaT for RAAN1 = ', RAAN1, ' and RAAN2 = ', RAAN2, ' is ', dt/86400, ' days.')
        return dt
    else:
        dt=-(2/3)*(RAAN1-RAAN2)/(constants.J2*(constants.R_EARTH**2)*((np.sqrt(constants.mu_EARTH/(a1**3))*np.cos(i1)/(a1**2))-(np.sqrt(constants.mu_EARTH/(a2**3))*np.cos(i2)/(a2**2))))
        if print_result == True:
            print('DeltaT for RAAN1 = ', RAAN1, ' and RAAN2 = ', RAAN2, ' is ', dt/86400, ' days.')
        return dt

def compute_dt_matrix(debris_data):

	""" 
    Function computing the delta_t matrix of all possible transfer from a debris to another

            inputs : 
            ------
                    - debris_data : dataframe
                            Orbital parameters of the debris considered
            output : 
            ------
                    - dt_matrix : array
                            Matrix whose (i,j) indice represents the delta_t required to modify the RAAN of the orbit from the RAAN of the debris i to the RAAN of the debris j

    """
	dt_matrix = np.zeros((constants.N_DEBRIS, constants.N_DEBRIS))
	for l in range(constants.N_DEBRIS):
		a=debris_data.values[l][0]
		i=debris_data.values[l][2]
		RAAN=debris_data.values[l][3]
		dt_matrix[l][l]=0.
		for j in range (l):
			dt=compute_dt(a,debris_data.values[j][0], i, debris_data.values[j][2], RAAN, debris_data.values[j][3])
			dt_matrix[l][j]=dt
			dt_matrix[j][l]=dt
	return dt_matrix


if __name__ == "__main__":

	# Setting the font

    rc('font', **{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)

	###############################
	#  							  #
	#  	    Recovering data       #
	#							  #
	###############################

    debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
    print(debris_data)

	###############################
	#  							  #
	#  	      Computing Δt        #
	#							  #
	###############################

    dt_matrix=compute_dt_matrix(debris_data)

	###############################
	#  							  #
	#    Plotting the results     #
	#							  #
	###############################

    dt_max=360
    ticks=np.arange(1,constants.N_DEBRIS+1,2)
    extent=(0.5,constants.N_DEBRIS+0.5,0.5,constants.N_DEBRIS+0.5)
    cmap0 = LinearSegmentedColormap.from_list('', ['darkblue', 'white'])
    plt.imshow(dt_matrix/86400,vmin=0, vmax=dt_max,interpolation='nearest',extent=extent, origin='lower', cmap=cmap0)
    plt.xlabel('Debris ID')
    plt.ylabel('Debris ID')
    plt.xticks(ticks)
    plt.yticks(ticks)
    colorbar=plt.colorbar()
    colorbar.set_label('$ \Delta t $ for RAAN alignment [days]')

    plt.show()
















