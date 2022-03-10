# -*- coding: utf-8 -*-
"""
Created on 16/12/2021

@author: Dev. Tom SEMBLANET
"""

import sys
import numpy as np
from utils import debris_data_loader as DDL

from regroupement.dV_computations.dV_matrix import dV_matrix_generation

from regroupement.optimizer.energy_computation import energy_computation	
from regroupement.optimizer.Recuit import Recuit

from utils.constants import COLORS


def script_generator(folder_path, coe_i_H, db_indexs_mat, debris_data_db):
	""" Generates the `Spacecraft` part of the 'MISSION_X' script (X = {1, 2, 3, 4, 5}, depending on the number
		of debris indexs passed in the `db_indexs` parameter). Once this script has been runned, the 'MISSION_X' 
		script is runnable in GMAT and simulate the trajectory of the Hunter S/C from debris to debris.

		inputs : 
        ------
            - folder_path: string
            	Path of the FOLDER in which the GMAT script is placed 
            - coe_i_H: float array
            	Initial orbital elements of the Hunter S/C
            - db_indexs_mat: list of list
            	Debris indexs for each groups
            - debris_data_db: <dataframe>
            	Orbital parameters of all debris

	
	"""

	for i, db_indexs in enumerate(db_indexs_mat):

		added_dt, total_dt = ['%----------------------------------------\n', '%---------- Spacecrafts\n', \
						'%----------------------------------------\n'], []

		added_dt.append("Create Spacecraft HUNTER;\n")
		added_dt.append("GMAT HUNTER.DateFormat = UTCGregorian;\n")
		added_dt.append("GMAT HUNTER.Epoch = '04 Jul 2014 09:30:15.235';\n")
		added_dt.append("GMAT HUNTER.CoordinateSystem = EarthMJ2000Eq;\n")
		added_dt.append("GMAT HUNTER.DisplayStateType = Keplerian;\n")
		added_dt.append("GMAT HUNTER.SMA = {};\n".format(coe_i_H[0]))
		added_dt.append("GMAT HUNTER.ECC = {};\n".format(coe_i_H[1]))
		added_dt.append("GMAT HUNTER.INC = {};\n".format(coe_i_H[2] * 180 / np.pi))
		added_dt.append("GMAT HUNTER.RAAN = {};\n".format(coe_i_H[3] * 180 / np.pi))
		added_dt.append("GMAT HUNTER.AOP = {};\n".format(coe_i_H[4] * 180 / np.pi))
		added_dt.append("GMAT HUNTER.TA = {};\n".format(coe_i_H[5] * 180 / np.pi))
		added_dt.append("GMAT HUNTER.OrbitColor = Red;\n")
		added_dt.append("\n")

		for k, indx in enumerate(db_indexs):
			added_dt.append("Create Spacecraft TARGET{};\n".format(k+1))
			added_dt.append("GMAT TARGET{}.DateFormat = UTCGregorian;\n".format(k+1))
			added_dt.append("GMAT TARGET{}.Epoch = '04 Jul 2014 09:30:15.235';\n".format(k+1))
			added_dt.append("GMAT TARGET{}.CoordinateSystem = EarthMJ2000Eq;\n".format(k+1))
			added_dt.append("GMAT TARGET{}.DisplayStateType = Keplerian;\n".format(k+1))
			added_dt.append("GMAT TARGET{}.SMA = {};\n".format(k+1, debris_data_db.values[indx, 0]))
			added_dt.append("GMAT TARGET{}.ECC = {};\n".format(k+1, debris_data_db.values[indx, 1]))
			added_dt.append("GMAT TARGET{}.INC = {};\n".format(k+1, debris_data_db.values[indx, 2] * 180 / np.pi))
			added_dt.append("GMAT TARGET{}.RAAN = {};\n".format(k+1, debris_data_db.values[indx, 3] * 180 / np.pi))
			added_dt.append("GMAT TARGET{}.AOP = {};\n".format(k+1, debris_data_db.values[indx, 4] * 180 / np.pi))
			added_dt.append("GMAT TARGET{}.TA = {};\n".format(k+1, debris_data_db.values[indx, 5] * 180 / np.pi))
			added_dt.append("GMAT TARGET{}.OrbitColor = {};\n".format(k+1, COLORS[k]))
			added_dt.append("\n")


		permanent_dt = []
		# Loading mission informations from model files
		with open(folder_path+"/models/MISSION_{}D_model".format(len(db_indexs)), "r") as file:
			dt = file.readlines()
			permanent_dt = dt

		total_dt = added_dt + permanent_dt
			
		with open(folder_path+"/MISSION_{}".format(i), "w") as file:
			for line in total_dt:
				file.write(line)

if __name__ == '__main__':

	folder_path = '/Users/semblanet/Desktop/Git/PIE-COS05/GMAT/MISSIONS/'

	coe_i_H = [7200, 0.0001, 0.9, 0.2, 0.4, 0.6]

	db_indexs_mat = [[25, 33, 24], [48, 8, 35], [29, 39, 11], [10, 35, 32], [30, 1, 16], [19, 6, 40], \
	[42, 43, 45], [41, 3, 23], [46, 0, 20], [47, 4, 49], [21, 5, 24], [26, 13, 31], [17, 22, 14], [44, 27, 36], \
	[28, 2, 38], [15, 7, 12]]

	# Loading of the debris database
	debris_data_db = DDL.convertTLEtoDF(DDL.recoveringDebrisData())

	script_generator(folder_path=folder_path, coe_i_H=coe_i_H, db_indexs_mat=db_indexs_mat, debris_data_db=debris_data_db)



