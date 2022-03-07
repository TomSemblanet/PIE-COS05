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


if __name__ == "__main__":

	# Loading debris data
	debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
	print('Debris data loaded')
	nb_debris = len(debris_data)

	# Implementing parameters for recuit/gibbs
	s_min = 4
	s_max = 5

	# Defining temperature
	# Ti = 1.5
	# Ti = 5
	# Tf = 0.001

	T = 0

	alpha = 0.95

	t_iter = 1000
	n_iter = 1

	n_classes = 500

	# Computing matrices
	DV = dV_matrix_generation(debris_data)
	DT = compute_dt_matrix(debris_data)

	E_min = sum(sum(DV + DT))

	V_tol = 1.0
	t_tol = 3*365.0

	n_try = 200

	# Launching Gibbs
	for i in range(n_try):
		print(i+1, ' / ', n_try)
		G_out, E_out, freq = Gibbs(nb_debris, DV, DT, T, s_min, s_max, n_classes, t_iter, n_iter, V_tol = V_tol, t_tol = t_tol)

		if E_out < E_min:
			E_min = E_out
			G_opt = G_out

	print('E_min = ', E_min)

	# E, E_transfers, E_transfers_dV, E_transfers_dt = energy_computation(G_opt, DV, DT, detail = True, cost = cost)
	# print(E_min)
	# print('detail of DV : ', E_transfers_dV)
	# print()
	# print('detail of Dt : ', E_transfers_dt)








