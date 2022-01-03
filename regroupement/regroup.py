import numpy as np
from numpy import linalg as LA

from . import recoveringDebrisData2 as RDB

from .utils import calcul_V

from .calculs_dV import SMA_dV, INC_dV, AOP_dV

debris_data = RDB.convertTLEtoDF(RDB.recoveringDebrisData())
ordered_debris = debris_data.sort_values(by=["RAAN (rad)"], ascending=False)


def deltaV_tot(i, j):
    a, e, i, omega, w, M, mass = ordered_debris.values[i][0], ordered_debris.values[i][1], ordered_debris.values[i][
        2], ordered_debris.values[i][3], ordered_debris.values[i][4], ordered_debris.values[i][5], ordered_debris.values[i][7]

    a2, e2, i2, omega2, w2, M2, mass2 = ordered_debris.values[j][0], ordered_debris.values[j][1], ordered_debris.values[j][
        2], ordered_debris.values[j][3], ordered_debris.values[j][4], ordered_debris.values[j][5], ordered_debris.values[j][7]

    # Coût pour la manoeuvre sur i :
    V1 = calcul_V(a)
    delta_V_i = INC_dV(i, i2, V1)
    # print(delta_V_i)
    i = i2

    # Coût pour la manoeuvre sur w :

    delta_V_w = AOP_dV(w, w2, omega, a, e, i, mass)
    # print(delta_V_w)
    w = w2

    # Coût pour la manoeuvre sur a :

    delta_V_a = SMA_dV(a, a2)
    # print(delta_V_a)
    a = a2

    delta_V_tot = delta_V_i + delta_V_w + delta_V_a
    return delta_V_tot


dict = {}
deltaV_matrix = np.zeros((len(ordered_debris), len(ordered_debris)))
for i in range(len(ordered_debris)):
    for j in range(i+1, len(ordered_debris)):
        dict[f"Total delta V between debris number {i} and {j} is "] = deltaV_tot(
            i, j)
        deltaV_matrix[i][j] = deltaV_tot(i, j)


print(dict)
