# -*- coding: utf-8 -*-
"""
Created on 16/12/2021

@author: Tom SEMBLANET
"""

import sys
import numpy as np
from utils import debris_data_loader as DDL

from regroupement.dV_computations.dV_matrix import dV_matrix_generation

from regroupement.optimizer.energy_computation import energy_computation	
from regroupement.optimizer.Recuit import Recuit

from utils.constants import COLORS


def script_generator(folder_path, coe_i_H, db_indexs):
	""" Generates the `Spacecraft` part of the 'MISSION_X' script (X = {1, 2, 3, 4, 5}, depending on the number
		of debris indexs passed in the `db_indexs` parameter). Once this script has been runned, the 'MISSION_X' 
		script is runnable in GMAT and simulate the trajectory of the Hunter S/C from debris to debris.

		inputs : 
        ------
            - folder_path: string
            	Path of the FOLDER in which the GMAT script is placed 
            - coe_i_H: float array
            	Initial orbital elements of the Hunter S/C
            - db_indexs: int array
            	Debris indexs 
	
	"""

	# Loading of the debris database
	debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())

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
		added_dt.append("GMAT TARGET{}.SMA = {};\n".format(k+1, debris_data.values[indx, 0]))
		added_dt.append("GMAT TARGET{}.ECC = {};\n".format(k+1, debris_data.values[indx, 1]))
		added_dt.append("GMAT TARGET{}.INC = {};\n".format(k+1, debris_data.values[indx, 2] * 180 / np.pi))
		added_dt.append("GMAT TARGET{}.RAAN = {};\n".format(k+1, debris_data.values[indx, 3] * 180 / np.pi))
		added_dt.append("GMAT TARGET{}.AOP = {};\n".format(k+1, debris_data.values[indx, 4] * 180 / np.pi))
		added_dt.append("GMAT TARGET{}.TA = {};\n".format(k+1, debris_data.values[indx, 5] * 180 / np.pi))
		added_dt.append("GMAT TARGET{}.OrbitColor = {};\n".format(k+1, COLORS[k]))
		added_dt.append("\n")


	permanent_dt = []
	# Deleting expired informations
	with open(folder_path+"/MISSION_{}".format(len(db_indexs)), "r") as file:
		start_indx = 0
		dt = file.readlines()

		for k, l in enumerate(dt):
			if '%---------- ForceModels' in l:
				start_indx = k
				break

		permanent_dt = dt[start_indx-1:]

	with open(folder_path+"/MISSION_{}".format(len(db_indexs)), "r") as file:
		total_dt = added_dt + permanent_dt
		
	with open(folder_path+"/MISSION_{}".format(len(db_indexs)), "w") as file:
		for line in total_dt:
			file.write(line)

if __name__ == '__main__':

	# Arguments to pass : GMAT SCRIPT FOLDER PATH + HUNTER INITIAL COE + DEBRIS INDEXS
	# -----------------
	script_generator(folder_path=sys.argv[1], coe_i_H=[float(arg) for arg in sys.argv[2:8]], db_indexs=[int(arg) for arg in sys.argv[8:]])



