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

	debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
	debris_data = debris_data.sort_values(by=["RAAN (rad)"], ascending=False)
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
	s_max = 5

	# Defining temperature
	# Ti = 1.5
	Ti = 2
	Tf = 0.001

	alpha = 0.97

	t_iter = 1000
	n_iter = 1
	# t_ier = 350
	# n_iter = 50

	#################################################################

	#################################################################
	#################################################################
	#							  CASE 0							#
	#################################################################
	#################################################################

	# case = 0 	# Case 'only_dV' i.e. we only consider the cost of the maneuvers

	#################################################################
	#################################################################
	#							  CASE 1							#
	#################################################################
	#################################################################

	case = 1	# Case 'dV_and_drift' i.e. we consider the cost of the maneuvers and the duration of the RAAN alignment

	#################################################################
	#################################################################
	#							  CASE 2							#
	#################################################################
	#################################################################

	# case = 2 	# Case 'only_dV' but taking in account maneuvers over the RAAN

	#################################################################


	if case == 0:

		n_classes = 10

		# Computing DV matrix
		DV = dV_matrix_generation(debris_data)
		DT = compute_dt_matrix(debris_data)

	elif case == 1:

		n_classes = 500

		# Computing DV and DT matrices
		DV = dV_matrix_generation(debris_data)
		DT = compute_dt_matrix(debris_data)

	else:

		n_classes = 100

		# Computing DV and DT matrices
		DV = dV_matrix_generation(debris_data, RAAN_maneuver = True)
		DT = compute_dt_matrix(debris_data)


	# Launching Recuit
	if case == 0 or case == 2:

		G_out, E_out, freqs = Recuit(nb_debris, s_min, s_max, DV, DT, Ti, Tf, alpha, n_classes, t_iter, n_iter)

		# Getting results
		E, E_transfers, E_transfers_dV = energy_computation(G_out, DV, DT, detail = True)

		# Expliciting detail of groups
		nb_grps = np.size(G_out,1)

		for k in range(nb_grps):
			grp = np.nonzero(G_out[:,k])[0]
			print('\nGroup number ', k+1 , ' = ', grp)
			print()

		print('\nFinal Energy = ', E_out, ' [km/s]')
		print('\n')
		print('Detail of delta_v = ', E_transfers_dV, ' [km/s]')
		print('\n')
		print('\nFinal state :\n', G_out)

	else:
		G_out, E_out, freqs = Recuit(nb_debris, s_min, s_max, DV, DT, Ti, Tf, alpha, n_classes, t_iter, n_iter, cost = 'dV_and_drift')

		# Getting results
		E, E_transfers, E_transfers_dV, E_transfers_dt = energy_computation(G_out, DV, DT, detail = True, cost = 'dV_and_drift')

		# Expliciting detail of groups
		nb_grps = np.size(G_out,1)

		for k in range(nb_grps):
			grp = np.nonzero(G_out[:,k])[0]
			print('\nGroup number ', k+1 , ' = ', grp)
			print()

		print('\nFinal Energy = ', E_out, ' [km/s]')
		print('\n')
		print('Detail of delta_v = ', E_transfers_dV, ' [km/s]')
		print('\n')
		print('Detail of elapsed time (J2) = ', E_transfers_dt, ' [days]')
		print('\nFinal state :\n', G_out)

		plot_bars(E_transfers_dV, np.array(E_transfers_dt)/365.0)






	# # launching Gibbs as decreasing gradient to evaluate local minimums
	# T = 0
	# G_out, E_out, freqs = Gibbs(nb_debris, DV, DT, T, s_min, s_max, n_classes, t_iter, n_iter)



