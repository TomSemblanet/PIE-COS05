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


def RAAN_evol(a,i,t,RAAN_0):
    """ 
    Function computing the delta_t [s] required to modify the RAAN of the orbit from RAAN1 to RAAN2 

    Arguments : 
        a (float) : Initial SMA [km]
        
        i (float) : Initial inclination [rad]

        t (float) : Time [s] at wich will be computed the RAAN from RAAN_0
                        
        RAAN_0 (float) : Initial Right ascension of the ascending node [rad]
        

    Returns : 
        RAAN_t (float) : Final Right ascension of the ascending node [rad]

    """
    n = np.sqrt(constants.mu_EARTH/(a**3))
    J2 = constants.J2
    R_eq = constants.R_EARTH

    RAAN_t = -3/2*n*((R_eq/a)**2)*J2*np.cos(i)*t + RAAN_0

    return RAAN_t


def compute_dt(a1,a2,i1,i2,RAAN1, RAAN2):
    """ 
    Function computing the delta_t [s] required to modify the RAAN of the orbit from RAAN1 to RAAN2 

    Arguments : 
        a1 (float) : Initial SMA [km]
        
        a2 (float) : Final SMA [km]
        
        i1 (float) : Initial inclination [rad]
        
        i2 (float) : Final inclination [rad]
        
        RAAN1 (float) : Initial Right ascension of the ascending node [rad]
        
        RAAN2 (float) : Final Right ascension of the ascending node [rad]
        

    Returns : 
        dt_days (float) : Required delta_t [days]

    """

    dt = 10
    RAAN_1_dt = RAAN_evol(a1,i1,dt,RAAN1)
    RAAN_2_dt = RAAN_evol(a2,i2,dt,RAAN2)

    delta_0 = np.abs(RAAN1 - RAAN2)
    delta_t = np.abs(RAAN_1_dt - RAAN_2_dt)

    if delta_0 > delta_t :
        delta_0 = np.abs(RAAN1 - RAAN2)
    else :
        delta_0 = np.abs(2*np.pi - delta_0)

    n1 = np.sqrt(constants.mu_EARTH/(a1**3))
    n2 = np.sqrt(constants.mu_EARTH/(a2**3))
    J2 = constants.J2
    R_eq = constants.R_EARTH

    denom = J2*(R_eq**2)*(n1*np.cos(i1)/(a1**2) - n2*np.cos(i2)/(a2**2))

    dt = np.abs(2.0/3*delta_0/denom)
    dt_days = dt/86400.0

    return dt_days

def compute_dt_matrix(debris_data):

    """ 
    Function computing the delta_t matrix of all possible transfer from a debris to another

    Arguments :         
        debris_data (dataframe) : Orbital parameters of the debris considered
        
    Returns :           
        dt_matrix (array) : Matrix whose (i,j) indice represents the delta_t (in seconds) required to modify the RAAN of the orbit from the RAAN of the debris i to the RAAN of the debris j

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
    #                             #
    #       Recovering data       #
    #                             #
    ###############################

    debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
    # print(debris_data)

    ###############################
    #                             #
    #         Computing Δt        #
    #                             #
    ###############################

    dt_matrix=compute_dt_matrix(debris_data)

    ###############################
    #                             #
    #    Plotting the results     #
    #                             #
    ###############################

    dt_max=360
    ticks=np.arange(1,constants.N_DEBRIS,2)
    extent=(0.5,constants.N_DEBRIS+0.5,0.5,constants.N_DEBRIS+0.5)
    cmap0 = LinearSegmentedColormap.from_list('', ['darkblue', 'white'])
    plt.imshow(dt_matrix,vmin=0, vmax=dt_max,interpolation='nearest',extent=extent, origin='lower', cmap=cmap0)
    plt.xlabel('Debris ID')
    plt.ylabel('Debris ID')
    plt.xticks(ticks)
    plt.yticks(ticks)
    colorbar=plt.colorbar()
    colorbar.set_label('$ \Delta t $ for RAAN alignment [days]')

    plt.show()














