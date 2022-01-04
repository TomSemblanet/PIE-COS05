# -*- coding: utf-8 -*-
"""
Created on 16/12/2021

@author: Yvan GARY
"""

from regroupement.delta_v_computation import recoveringDebrisData2 as RDB
from regroupement.delta_v_computation.DV_computation import Generate_DV_Matrix

from regroupement.optimizer.Recuit import Recuit


if __name__ == "__main__":

	debris_data = RDB.convertTLEtoDF(RDB.recoveringDebrisData())

	# Implementing arguments for Recuit
	nb_debris = 19
	card_grp = 5

	DV = Generate_DV_Matrix(debris_data)

	Ti = 0.1
	Tf = 0.01

	alpha = 0.95

	n_classes = 10
	t_iter = 100
	n_iter = 50

	# Launching Recuit
	G_out, E_out, freqs = Recuit(nb_debris, card_grp, DV, Ti, Tf, alpha, n_classes, t_iter, n_iter)

	print('Final global delta_v = ', E_out)
	print('\n')
	print('Final state :\n', G_out)