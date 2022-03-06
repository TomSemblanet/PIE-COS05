# -*- coding: utf-8 -*-
"""
Created on 09/12/2021

@author: Yvan GARY
"""

import numpy as np
from itertools import permutations

def single_energy_computation(debris, DV, DT):
	''' Function used to compute the delta t (J2) associated to a state

	Arguments:
		debris (array) : Array containing a combination of debris 
		
		DV (Matrix): Matrix containing the delta_v associated to each maneuver		

		DT (Matrix): Matrix containing the delta_t associated to each J2 drift

	Returns: 
		dV (float): Global delta v associated to the state G
		
		dT (float): Global elapsed time for RAAN alignement

	'''

	dV = 0
	dT = 0
	nb_debris = len(debris)

	for l in range(nb_debris-1):
		dv = DV[debris[l],debris[l+1]]
		dt = DT[debris[l],debris[l+1]]

		dV += dv
		dT += dt

	return dV,dT

def minimal_energy_computation(ordered_debris, DV, DT, V_tol, t_tol):
	''' Function used to compute the energy associated to a state

	Arguments:
		ordered_debris (array) : array containing the debris in their natural order
		
		DV (Matrix): Matrix containing the delta_v associated to each maneuver
		
		DT (Matrix): Matrix containing the elapsed time associated to each "J2 perturbation duration" between two debris

		V_tol (float) : Order of magnitude of tolerated dV for one mission in km/s

		t_tol (float) : Order of magnitude of tolerated duration for a mission in days

	Returns: 
		E_min (float): Minimal energy produced by the optimal combination of debris 

		debris_opt (array) : Optimal debirs configuration for successive maneuvers
	'''

	E_min = 10**5
	debris_opt = ordered_debris

	for debris in permutations(ordered_debris):
		dV, dT = single_energy_computation(debris,DV,DT)

		# E = dV/V_tol + dT/t_tol
		E = dV + V_tol/t_tol*dT

		if E<E_min:
			E_min = E
			dV_opt = dV
			dT_opt = dT
			debris_opt = debris

	# return E_min, debris_opt
	return dV_opt, dT_opt, debris_opt



def energy_computation(G, DV, DT, V_tol, t_tol, show_grps = False):
	''' Function used to compute the energy associated to a state

	Arguments:
		G (Matrix): Actual state for which we compute the energy
		
		DV (Matrix): Matrix containing the delta_v associated to each maneuver
		
		DT (Matrix): Matrix containing the elapsed time associated to each "J2 perturbation duration" between two debris

		V_tol (float) : Order of magnitude of tolerated dV for one mission in km/s

		t_tol (float) : Order of magnitude of tolerated duration for a mission in days
		
		show_grps (bool) : False by default - If True, gives the detail of groups formed in the state

	Returns: 
		E (float): Energy associated to the state G
		
		grps (array) - optionnal : False by default. If true, gives the detail of the order of the maneuvers for each group

	'''

	nb_grp = np.size(G,1)
	if show_grps:
		grps = []

	E = 0

	for n in range(nb_grp):

		grp = G[:,n]
		ordered_debris = np.nonzero(grp)[0]

		# E_min, debris_opt = minimal_energy_computation(ordered_debris,DV,DT,V_tol,t_tol)
		dV, dT, debris_opt = minimal_energy_computation(ordered_debris,DV,DT,V_tol,t_tol)

		std = np.var([dV/V_tol, dT/t_tol])
		E_n = dV/V_tol + dT/t_tol + std

		E += E_n
		if show_grps:
			grps.append(debris_opt)

	if show_grps:
		return E, grps
	else:
		return E


