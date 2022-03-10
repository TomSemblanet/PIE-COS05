import numpy as np
import matplotlib.pyplot as plt

from utils import debris_data_loader as DDL
import regroupement.dV_computations.compute_dt_alignment as cdt
from regroupement.dV_computations.dV_matrix import dV_matrix_generation
from regroupement.optimizer.energy_computation import single_energy_computation


if __name__ == "__main__":

	print('Loading debris data...')
	debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
	print(debris_data)

	DV = dV_matrix_generation(debris_data)

	debris = [44,39,11]
	# a1,a2 = debris_data.values[debris,0]
	# i1,i2 = debris_data.values[debris,2]
	# raan1,raan2 = debris_data.values[debris,3] 


	# dt1 = cdt.compute_dt(a1,a2,i1,i2,raan1,raan2)
	# dt2 = cdt.compute_dt(a1,a2,i1+0.3*np.pi/180,i2+0.3*np.pi/180,raan1,raan2)

	# print(a1)
	# print(a2)
	# print(i1*180/(np.pi))
	# print(i2*180/(np.pi))
	# print(raan1*180/(np.pi))
	# print(raan2*180/(np.pi))

	# print()

	# print(dt1)
	# print(dt2)

	dv,dT,dt_ = single_energy_computation(debris,DV,debris_data)

	print('Deuxième liste: ', dt_)