import numpy as np

from utils import debris_data_loader as DDL
from utils.constants import mu_EARTH

from regroupement.dV_computations.maneuvers_dV import SMA_dV, INC_dV, AOP_dV


def total_dV_computation(ordered_debris, k, l):
    """ 
    Function computing the delta_v value [km/s] between the orbits of the debris k and l.

    Arguments:
        ordered_debris (DataFrame): Data of the debris considered sorted by decreasing w
        k (int) : label of the first debris
        l (int) : label of the second debris

    Returns:
        delta_V_tot (float) : delta_v corresponding to the transfer from orbit of debris k to orbit of debris l [km/s]

    """

    # The RAAN of the target orbit has to be greater or equal to the RAAN of the initial orbit
    # Otherwize, the use of J2 effect isn't usable and we suppose that the transfer is infeasible
    if k > l :
        dV_TOT = 0

    else:
        # Extraction of the initial and final orbits orbital elements (+ mass)
        a, e, i, omega, w, M, mass = ordered_debris.values[k][[0, 1, 2, 3, 4, 5, 8]]
        a2, e2, i2, omega2, w2, M2, mass2 = ordered_debris.values[l][[0, 1, 2, 3, 4, 5, 8]]

        # Inclination modification cost [m/s]
        V1 = np.sqrt(mu_EARTH/a)
        dV_INC = INC_dV(i, i2, V1)
        i = i2

        # AOP modification cost [m/s]
        dV_AOP = AOP_dV(w, w2, omega, a, e, i, mass)
        w = w2

        print("{} km/s".format(dV_AOP))

        # SMA modification cost [m/s]
        dV_SMA = SMA_dV(a, a2)
        a = a2

        # Total delta-V [m/s]
        dV_TOT = dV_INC + dV_AOP + dV_SMA

    return dV_TOT

def dV_matrix_generation(debris_data):
    """ 
    Function computing the delta_v matrix of all possible transfer from a debris to another

    Arguments:
        debris_data (DataFrame): Data of the debris considered

    Returns:
        dV_matrix (Matrix) : delta-V matrix (triangular superior) recording all possible values of delta_v for the given set of debris
    """

    # Ordering debris according to decreasing RAAN
    ordered_debris = debris_data.sort_values(by=["RAAN (rad)"], ascending=False)

    dV_matrix = np.zeros((len(ordered_debris), len(ordered_debris)))

    for i in range(len(ordered_debris)):
         for j in range(i+1, len(ordered_debris)):
            dV_matrix[i][j] = total_dV_computation(ordered_debris, i, j)

    return dV_matrix

if __name__ == '__main__':
    
    # Loading and sorting of the debris data
    debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
    
    # Computation of the delta-V matrix 
    dV_matrix = dV_matrix_generation(debris_data)
