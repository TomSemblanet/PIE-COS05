# -*- coding: utf-8 -*-
"""
Created on 16/12/2021

@author: Yvan GARY
"""

import numpy as np

from utils import debris_data_loader as DDL

from regroupement.dV_computations.dV_matrix import dV_matrix_generation
from regroupement.dV_computations.compute_dt_alignment import compute_dt_matrix

from regroupement.optimizer.energy_computation import energy_computation	
from regroupement.optimizer.Recuit import Recuit
from regroupement.optimizer.Gibbs import Gibbs


if __name__ == "__main__":

	debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
	print(debris_data)

	# Implementing arguments for Recuit
	nb_debris = len(debris_data)

	# Max and min number of debris per group
	s_min = 4
	s_max = 5

	# Computing DV and DT matrices
	DV = dV_matrix_generation(debris_data)
	DT = compute_dt_matrix(debris_data)

	# Ti = 0.5
	# Tf = 0.001

	Ti = 1
	Tf = 0.01

	T=0

	alpha = 0.95

	n_classes = 100
	t_iter = 150
	n_iter = 50

	# Launching Recuit
	G_out, E_out, freqs = Recuit(nb_debris, s_min, s_max, DV, DT, Ti, Tf, alpha, n_classes, t_iter, n_iter)

	# Getting results
	E, E_transfers, E_transfers_dV, E_transfers_dt = energy_computation(G_out, DV, DT, detail = True)

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


