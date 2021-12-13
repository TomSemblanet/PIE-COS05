#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 10:19:33 2021

@author: g.pierre
"""
# pip install -U TLE-tools
# pip install spacetrack
import spacetrack
from spacetrack import SpaceTrackClient
from tletools import TLE
import numpy as np
import pandas as pd


def recoveringDebrisData(*args):

    st = SpaceTrackClient('pierre.gaetan@outlook.com', 'zywraj-Tadky0-fezgek')
    if len(args) != 0:
        TLE_String = st.tle_latest(
            norad_cat_id=[args], ordinal=1, format='tle')
    else:
        TLE_String = st.tle_latest(norad_cat_id=[27001, 27601, 15334, 10732, 24279, 21090, 15772, 10693, 27387,
                                                 7594, 23180, 10138, 13917, 13719, 14625, 12092, 9044, 12504, 16292], ordinal=1, format='tle')
    return TLE_String


def convertTLEtoDF(TLE_String):

    columns = ['a (km)', 'e', 'i (rad)', 'RAAN (rad)', 'Argument of periapsis (rad)',
               'Mean anomaly (rad)', 'Date (year)', 'Date (fraction of year)']
    index = [27001, 27601, 15334, 10732, 24279, 21090, 15772, 10693, 27387,
             7594, 23180, 10138, 13917, 13719, 14625, 12092, 9044, 12504, 16292]
    TLE_DF = pd.DataFrame(index=index, columns=columns)
    mu = 398600.4418
    TLE_String = recoveringDebrisData()
    tle_lines = TLE_String.strip().splitlines()
    Current_TLE_String = []
    for i in range(int(len(tle_lines)/2)):
        Current_TLE_String = [' \n', tle_lines[2*i], tle_lines[2*i+1]]
        tle = TLE.from_lines(*Current_TLE_String)
        TLE_DF.at[index[i], 'a (km)'] = (
            mu/(((2*np.pi/86400)*tle.n)**2))**(1/3)
        TLE_DF.at[index[i], 'e'] = tle.ecc
        TLE_DF.at[index[i], 'i (rad)'] = np.radians(tle.inc)
        TLE_DF.at[index[i], 'RAAN (rad)'] = np.radians(tle.raan)
        TLE_DF.at[index[i], 'Argument of periapsis (rad)'] = np.radians(
            tle.argp)
        TLE_DF.at[index[i], 'Mean anomaly (rad)'] = np.radians(tle.M)
        TLE_DF.at[index[i], 'Date (year)'] = tle.epoch_year
        TLE_DF.at[index[i], 'Date (fraction of year)'] = tle.epoch_day

    return TLE_DF
