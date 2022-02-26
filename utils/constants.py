import numpy as np 
import math as mt

mu_EARTH = 398600.4418 # Earth's gravitational parameter [km^3/s^-2]
R_EARTH = 6384.415 # Earth's radius [km]
J2 = 0.001082 # J2 Perturbation

NORAD_ID_DEBRIS_LESS_3T=[27001, 27601, 15334, 10732, 24279, 21090, 15772, 10693, 27387, \
 7594, 23180, 10138, 13917, 13719, 14625, 12092, 9044, 12504, 16292] # NORAD_ID of the debris we chose

NORAD_ID_DEBRIS=[22566, 22220, 31793, 26070, 16182, 20625, 27006, 23705, 25407, 23405, 17974,\
23088, 22285,22803, 19650, 24298, 28353, 17590, 19120, 25400, 27386, 27001, 24277, 27601, 15334, \
37932, 10732, 24279, 23704, 21090, 28352, 23087, 19119, 27597, 25861, 15772, 10693, 17973, 27387, \
7594, 23180, 10138, 13917, 13719, 14625, 20624, 12092, 9044, 12504, 16292]

MASSES_DEBRIS=[9000, 9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,\
9000,9000,9000,9000,9000,7800, 2500, 3560, 3000, 2440, 4000, 1435, 3250, 1435, 3250, 3250, 3250, 3680, \
9000, 2440, 1435, 3250, 2575, 1435, 1435, 1435, 1435, 1100, 1435, 3250, 1435, 1435, 800, 1435]

N_DEBRIS=len(NORAD_ID_DEBRIS) # Number of debris