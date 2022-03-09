# -*- coding: utf-8 -*-
"""
Created on 16/12/2021

@author: Yvan GARY
"""

import numpy as np
import matplotlib.pyplot as plt

from utils import debris_data_loader as DDL

from regroupement.dV_computations.dV_matrix import dV_matrix_generation
from regroupement.dV_computations.compute_dt_alignment import compute_dt_matrix

from regroupement.optimizer.energy_computation import energy_computation
from regroupement.optimizer.energy_computation import single_energy_computation

from regroupement.optimizer.Recuit import Recuit
from regroupement.optimizer.Gibbs import Gibbs

def plot_bars(Delta_v, Delta_t):

	nb_groups = len(Delta_v)

	debris = range(1,nb_groups+1)
	debris = np.array(debris)

	plt.bar(debris-0.2, Delta_v, color = 'r', width = 0.4, label = 'Delta V [km/s]')
	plt.bar(debris+0.2, Delta_t, color = 'g', width = 0.4, label = 'J2 Drift [years]')
	plt.legend()
	plt.xlabel('Groups')
	plt.show()


if __name__ == "__main__":

	print('Loading debris data...')
	debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
	print(debris_data)

	############### Implementing arguments for Recuit ############### 
	nb_debris = len(debris_data)

	################ PLOTTING OF ORBITAL PARAMETERS OF DEBRIS #################
	# ordered_derbis = debris_data.sort_values(by=["RAAN (rad)"], ascending=True)

	# RAANs = []
	# INCs = []
	# SEMs = []

	# for k in range(nb_debris):
	# 	sem = debris_data.values[k][0]
	# 	raan = debris_data.values[k][3]
	# 	inc = debris_data.values[k][2]
	# 	SEMs.append(sem)
	# 	RAANs.append(raan)
	# 	INCs.append(inc)

	# # plt.plot(range(nb_debris), RAANs, '+', color = 'green')
	# # plt.plot(range(nb_debris), INCs, '+', color = 'red')
	# plt.plot(range(nb_debris), SEMs, '+', color = 'red')
	# plt.xlabel('ordered debris')
	# plt.ylabel('SEM')
	# # plt.title('Increasing values of RAAN and INC values')
	# plt.show()

	# input()

	###########################################################################

	# Max and min number of debris per group
	s_min = 4
	s_max = 4

	# Defining temperature
	Ti = 0.001
	Tf = 0.00098

	alpha = 0.95

	t_iter = 500
	n_iter = 1
	# t_ier = 350
	# n_iter = 50

	n_classes = 150

	# Computing DV and DT matrices
	# DV = dV_matrix_generation(debris_data)
	# DT = compute_dt_matrix(debris_data)

	# Tolerances
	# V_tol = 1.0
	V_tol = 1.0
	t_tol = 365.0

	# CHOIX DU SCENARIO
	scenario = 1
	print()
	print('Scenario = ', scenario)

	# Scenario 1
	if scenario == 1:
		DV = dV_matrix_generation(debris_data)
		G_out, E_out, freqs = Recuit(nb_debris, s_min, s_max, DV, debris_data, Ti, Tf, alpha, n_classes, t_iter, n_iter, V_tol = V_tol, t_tol = t_tol)

		E, grps = energy_computation(G_out, DV, debris_data, V_tol = V_tol, t_tol = t_tol, show_grps = True)

	# Scénario 2
	elif scenario == 2:
		DV = dV_matrix_generation(debris_data, RAAN_maneuver = True)
		G_out, E_out, freqs = Recuit(nb_debris, s_min, s_max, DV, debris_data, Ti, Tf, alpha, n_classes, t_iter, n_iter, V_tol = 0.0, t_tol = t_tol)

		E, grps = energy_computation(G_out, DV, DT, V_tol = V_tol, t_tol = t_tol, show_grps = True)

	# Scénario 3
	elif scenario == 3:
		DV = dV_matrix_generation(debris_data)
		dV_matrix_RAAN = dV_matrix_generation(debris_data, RAAN_maneuver = True)
		G_out, E_out, freqs = Recuit(nb_debris, s_min, s_max, DV, debris_data, Ti, Tf, alpha, n_classes, t_iter, n_iter, V_tol = V_tol, t_tol = t_tol, DV_RAAN = dV_matrix_RAAN)

		E, grps, RAAN_mans = energy_computation(G_out, DV, debris_data, V_tol = V_tol, t_tol = t_tol, DV_RAAN = dV_matrix_RAAN, show_grps = True)

	else:
		raise Exception('This scenario has not been implemented, please choose between 1,2 or 3.')

	all_dV = []
	all_dT = []

	count = 0
	print('Groups : \n')
	for grp in grps:
		print('Group ', count + 1 , ' : ', grp)
		if scenario == 1 or scenario == 2:
			dV,dT = single_energy_computation(grp,DV,debris_data)
		else:
			dV,dT, RAAN_maneuver = single_energy_computation(grp,DV,debris_data,DV_RAAN=dV_matrix_RAAN, V_tol = V_tol, t_tol = t_tol)
			print('Associated RAAN maneuver : ', RAAN_mans[count])

		print()
		all_dV.append(dV)
		all_dT.append(dT)
		count += 1
		
	all_dV = np.array(all_dV)
	all_dT = np.array(all_dT)

	print('Mean Delta_v : ', np.mean(all_dV), '[km/s]')
	print()
	print('detail of delta V : ', all_dV)
	print()
	if scenario == 1 or scenario == 3 :
		print('Mean mission duration : ', np.mean(all_dT)/365.0, '[years]')
		print()
		print('detail of mission duration : ', all_dT)
		plot_bars(all_dV, all_dT/365.0)
	else:
		nb_groups = len(all_dV)

		debris = range(1,nb_groups+1)
		debris = np.array(debris)

		plt.bar(debris, all_dV, color = 'r', width = 0.4, label = 'Delta V [km/s]')
		plt.legend()
		plt.xlabel('Groups')
		plt.show()




