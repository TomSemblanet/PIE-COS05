# -*- coding: utf-8 -*-
"""
Created on 16/12/2021

@author: Yvan GARY
"""

from utils import debris_data_loader as DDL
from regroupement.dV_computations.dV_matrix import dV_matrix_generation

from regroupement.optimizer.Recuit import Recuit


if __name__ == "__main__":

	debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())

	# Implementing arguments for Recuit
	nb_debris = len(debris_data)
	card_grp = 5

	DV = dV_matrix_generation(debris_data)

	Ti = 1
	Tf = 0.005

	alpha = 0.95

	n_classes = 100
	t_iter = 100
	n_iter = 50

	# Launching Recuit
	G_out, E_out, freqs = Recuit(nb_debris, card_grp, DV, Ti, Tf, alpha, n_classes, t_iter, n_iter)

	print('Final global delta_v = ', E_out)
	print('\n')
	print('Final state :\n', G_out)