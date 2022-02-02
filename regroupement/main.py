# -*- coding: utf-8 -*-
"""
Created on 16/12/2021

@author: Yvan GARY
"""

from utils import debris_data_loader as DDL

from regroupement.dV_computations.dV_matrix import dV_matrix_generation

from regroupement.optimizer.energy_computation import energy_computation	
from regroupement.optimizer.Recuit import Recuit
from regroupement.optimizer.Gibbs import Gibbs


if __name__ == "__main__":

	debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
	print(debris_data)

	# Implementing arguments for Recuit
	nb_debris = len(debris_data)
	card_grp = 5

	# Max and min number of debris per group
	s_min = 4
	s_max = 5

	DV = dV_matrix_generation(debris_data)

	Ti = 0.5
	Tf = 0.001
	T=0

	alpha = 0.95

	n_classes = 10
	t_iter = 150
	n_iter = 50

	# Launching Recuit
	G_out, E_out, freqs = Recuit(nb_debris, s_min, s_max, DV, Ti, Tf, alpha, n_classes, t_iter, n_iter)

	# G_out, E_out, freqs = Gibbs(nb_debris, card_grp, DV, T, n_classes, t_iter, n_iter)

	E, E_transfers = energy_computation(G_out, DV, detail = True)

	print('\nFinal global delta_v = ', E_out, ' [km/s]')
	print('\n')
	print('Detail of delta_v = ', E_transfers, ' [km/s]')
	print('\nFinal state :\n', G_out)


