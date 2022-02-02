""" This module computes the delta-Vs necessary to change orbital parameters : 
        - SMA (semi-major axis) 
        - ECC (eccentricity)
        - INC (inclination)
        - AOP (argument of perigee)
        - RAAN (right ascension of the ascending node)
        - TA (true anomaly)
"""

import numpy as np

from utils.coc import kep2cart
from utils.constants import mu_EARTH


def SMA_dV(a1, a2):
    """ Computes the delta-V [km/s] required to modify the SMA of the orbit from ai to af 

            inputs : 
            ------
                    - a1 : float
                            Initial SMA [km]
                    - a2 : float
                            Final SMA [km]
            ouputs : 
            ------
                    - dV : float
                            Required delta-V [km/s]
    """

    # Speed on the orbit assuming it to be circular [km/s]
    V = np.sqrt(mu_EARTH/a1)
    delta_a = np.abs(a2 - a1)

    dV = V * delta_a / (2*a1)

    return dV


def INC_dV(i1, i2, V1):
    """ Computes the delta-V [km/s] required to modify the INC of the orbit from i1 to i2 

            inputs : 
            ------
                    - a1 : float
                            Initial INC [rad]
                    - a2 : float
                            Final INC [rad]
                    - V1 : float
                            Velocity on the orbit (supposed to be circular)
            ouputs : 
            ------
                    - dV : float
                            Required delta-V [km/s]
    """

    delta_i = np.abs(i2 - i1)
    dV = 2 * V1 * np.sin(delta_i / 2)

    return dV


def AOP_dV(w1, w2, RAAN, a, e, i, m):
    """ Computes the delta-V [km/s] required to modify the AOP of the orbit from w1 to w2 
            inputs : 
            ------
                    - w1 : float
                            Initial AOP [rad]
                    - w2 : float
                            Final AOP [rad]
                    - RAAN : float
                            Right ascension of the ascending node [rad]
                    - a : float
                            Semi-major axis [km]
                    - e : float
                            Eccentricity [-]
                    - i : float
                            Inclination [rad]
                    - m : float
                            Body's mass [kg]
            ouputs : 
            ------
                    - dV : float
                            Required delta-V [km/s]
    """

    delta_w = np.abs(w2 - w1)
    dV = 2 * e * np.sqrt(mu_EARTH / (a*(1 - e**2))) * np.sin(delta_w/2)

    return dV



