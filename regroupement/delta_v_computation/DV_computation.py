
import numpy as np
from numpy import linalg as LA

# from . import recoveringDebrisData2 as RDB
from regroupement.delta_v_computation import recoveringDebrisData2 as RDB

from regroupement.delta_v_computation.utils import calcul_V

from regroupement.delta_v_computation.calculs_dV import SMA_dV, INC_dV, AOP_dV

# debris_data = RDB.convertTLEtoDF(RDB.recoveringDebrisData())
# ordered_debris = debris_data.sort_values(by=["RAAN (rad)"], ascending=False)


def deltaV_tot(ordered_debris, k, l):

    ''' Function computing the delta_v value [km/s] between the orbits of the debris k and l.

    Arguments:
        ordered_debris (DataFrame): Data of the debris considered sorted by decreasing w
        k (int) : label of the first debris
        l (int) : label of the second debris

    Returns:
        delta_V_tot (float) : delta_v corresponding to the transfer from orbit of debris k to orbit of debris l [km/s]

    '''

    if k > l :
        # We consider no transfer from an orbit with lower w to higher w
        delta_V_tot = 0

    else :

        a, e, i, omega, w, M, mass = ordered_debris.values[k][0], ordered_debris.values[k][1], ordered_debris.values[k][
            2], ordered_debris.values[k][3], ordered_debris.values[k][4], ordered_debris.values[k][5], ordered_debris.values[k][7]

        a2, e2, i2, omega2, w2, M2, mass2 = ordered_debris.values[l][0], ordered_debris.values[l][1], ordered_debris.values[l][
            2], ordered_debris.values[l][3], ordered_debris.values[l][4], ordered_debris.values[l][5], ordered_debris.values[l][7]

        # Cost of the maneuver on i
        V1 = calcul_V(a)
        delta_V_i = INC_dV(i, i2, V1)
        i = i2

        # Cost of the maneuver on w

        delta_V_w = AOP_dV(w, w2, omega, a, e, i, mass)
        w = w2

        # Cost of the maneuver on a :

        delta_V_a = SMA_dV(a, a2)
        a = a2

        delta_V_tot = delta_V_i + delta_V_w + delta_V_a

    return delta_V_tot

def Generate_DV_Matrix(debris_data):

    ''' Function computing the delta_v matrix of all possible transfer from a debris to another

    Arguments:
        debris_data (DataFrame): Data of the debris considered

    Returns:
        deltaV_matrix (Matrix) : delta_v matrix (triangular superior) recording all possible values of delta_v for the given set of debris

    '''
    # Ordering debris according to decreasing w
    ordered_debris = debris_data.sort_values(by=["RAAN (rad)"], ascending=False)

    deltaV_matrix = np.zeros((len(ordered_debris), len(ordered_debris)))

    for i in range(len(ordered_debris)):
         for j in range(i+1, len(ordered_debris)):
            deltaV_matrix[i][j] = deltaV_tot(ordered_debris, i, j)

    return deltaV_matrix

