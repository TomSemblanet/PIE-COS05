import numpy as np
import matplotlib.pyplot as plt

from utils import debris_data_loader as DDL

from regroupement.dV_computations.dV_matrix import dV_matrix_generation
from regroupement.dV_computations.compute_dt_alignment import compute_dt_matrix

from regroupement.optimizer.energy_computation import energy_computation
from regroupement.optimizer.energy_computation import single_energy_computation
from regroupement.optimizer.energy_computation import minimal_energy_computation

if __name__ == "__main__":
	print('Loading debris data...')
	debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
	print('Debris data loaded')


	DV = dV_matrix_generation(debris_data)
	DV_RAAN = dV_matrix_generation(debris_data, RAAN_maneuver = True)
	DT = compute_dt_matrix(debris_data)

	grp = [13,18,0,44]
	RAAN_man = [0,0,0]

	N = len(grp)

	V_tol = 1.0
	t_tol = 3*365.0

	dV = 0
	dT = 0

	for n in range(N-1):
		if RAAN_man[n]:
			dv = DV_RAAN[grp[n],grp[n+1]]
			dV+=dv
		else:
			dv = DV[grp[n],grp[n+1]]
			dt = DT[grp[n],grp[n+1]]
			dV+=dv
			dT+=dt

	print(dV)
	print(dT)

	# dV_,dT_ = single_energy_computation(grp,DV,DT)

	# print('dV method single = ', dV_)
	# print('dT method single = ', dT_)






