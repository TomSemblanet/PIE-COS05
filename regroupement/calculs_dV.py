""" This module computes the delta-Vs necessary to change orbital parameters : SMA, ECC, INC, AOP, RAAN and TA """

import numpy as np 

from utilitaires.coc import kep2cart
from utilitaires.constants import mu_EARTH


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

	V = np.sqrt(mu/ai) # speed in the orbit assuming it to be circular [km/s]
	delta_a = np.abs(a2 - a1)

	dV = V * delta_a / 2 / a1

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

	dV = 2 * V1 * np.sin(np.abs(i2 - i1) / 2)

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
	ta = 0.5 * delta_w

	coe = np.array([a, e, i, w1, RAAN, ta]) 
	r = kep2cart(coe, mu_EARTH)

	h = np.linalg.norm(np.cross(r[:3], r[3:]))

	dV = 2 * mu_EARTH * m * e * np.sin(delta_w/2) / L

	return dV
