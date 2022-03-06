#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 10:19:33 2021

@author: g.pierre
"""
## pip install -U TLE-tools
## pip install spacetrack 

from tletools import TLE
import numpy as np
import pandas as pd
import spacetrack
from utils import constants
from spacetrack import SpaceTrackClient

def recoveringDebrisData(*args):   
    """ Function which aims to recover data from orbiting objects from Space-Track.org 
    
    Arguments : 
    
        args (list of int) : Norad IDs of objects of interest, default value is the list given in constant.py

    Returns : 
    
        TLE_String (string) : Concatenated TLEs (string) of each object
    """
    
    st = SpaceTrackClient('pierre.gaetan@outlook.com', 'zywraj-Tadky0-fezgek')
    if len(args)!=0:
        TLE_String=st.tle_latest(norad_cat_id=[args], ordinal=1, format='tle')                   
    else:
        TLE_String=st.tle_latest(norad_cat_id=constants.NORAD_ID_DEBRIS, ordinal=1, format='tle')           
    return TLE_String

def convertTLEtoDF(TLE_String):    
    """ Function which aims to convert data in TLE format towards pandas dataframe (for our set of debris) 
    
    Arguments : 
    
        TLE_String (string) : Concatenated TLEs (string) of each object

    Returns : 
    
        TLE_DF (pandas dataframe) : DataFrame containing orbital parameters, time and mass for each debris
    """
    
    columns=[ 'a (km)', 'e', 'i (rad)', 'RAAN (rad)', 'Argument of periapsis (rad)', 'Mean anomaly (rad)', 'Date (year)', 'Date (fraction of year)', 'Mass (kg)']
    index=constants.NORAD_ID_DEBRIS
    TLE_DF=pd.DataFrame(index=index, columns=columns);
    mass=constants.MASSES_DEBRIS
    tle_lines = TLE_String.strip().splitlines()
    Current_TLE_String=[]
    for i in range(int(len(tle_lines)/2)):
        Current_TLE_String=[' \n', tle_lines[2*i], tle_lines[2*i+1]]
        tle = TLE.from_lines(*Current_TLE_String)
        TLE_DF.at[index[i], 'a (km)']= (constants.mu_EARTH/(((2*np.pi/86400)*tle.n)**2))**(1/3)
        TLE_DF.at[index[i], 'e']= tle.ecc
        TLE_DF.at[index[i], 'i (rad)']= np.radians(tle.inc)
        TLE_DF.at[index[i], 'RAAN (rad)']= np.radians(tle.raan)
        TLE_DF.at[index[i], 'Argument of periapsis (rad)']= np.radians(tle.argp)
        TLE_DF.at[index[i], 'Mean anomaly (rad)']= np.radians(tle.M)
        TLE_DF.at[index[i], 'Date (year)']= tle.epoch_year
        TLE_DF.at[index[i], 'Date (fraction of year)']= tle.epoch_day
        TLE_DF.at[index[i], 'Mass (kg)']= mass[i]
        
    return TLE_DF 
    
