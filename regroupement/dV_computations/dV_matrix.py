import numpy as np

from utils import debris_data_loader as DDL
from utils.constants import mu_EARTH
from utils import constants
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from regroupement.dV_computations.maneuvers_dV import SMA_dV, INC_dV, AOP_dV, RAAN_dV


def total_dV_computation(debris, k, l, RAAN_maneuver = False):
    """ 
    Function computing the delta_v value [km/s] between the orbits of the debris k and l.

    Arguments:
        debris (DataFrame): Data of the debris considered sorted by decreasing w
        
        k (int) : label of the first debris
        
        l (int) : label of the second debris

        RAAN_maneuver (boolean) - False by default : Boolean indicating if RAAN maneuvers are taken in account in the delta V computation
        

    Returns:
        delta_V_tot (float) : delta_v corresponding to the transfer from orbit of debris k to orbit of debris l [km/s]

    """

    # The RAAN of the target orbit has to be greater or equal to the RAAN of the initial orbit
    # Otherwize, the use of J2 effect isn't usable and we suppose that the transfer is infeasible
    if k > l :
        dV_TOT = 0

    else:
        # Extraction of the initial and final orbits orbital elements (+ mass)
        a, e, i, omega, w, M, mass = debris.values[k][[0, 1, 2, 3, 4, 5, 8]]
        a2, e2, i2, omega2, w2, M2, mass2 = debris.values[l][[0, 1, 2, 3, 4, 5, 8]]

        # Inclination modification cost [m/s]
        V1 = np.sqrt(mu_EARTH/a)
        dV_INC = INC_dV(i, i2, V1)
        i = i2

        # AOP modification cost [m/s]
        dV_AOP = AOP_dV(w, w2, omega, a, e, i, mass)
        w = w2

        #print("{} km/s".format(dV_AOP))

        # SMA modification cost [m/s]
        dV_SMA = SMA_dV(a, a2)
        a = a2

        # Total delta-V [m/s]
        if RAAN_maneuver:
            dV_RAAN = RAAN_dV(omega, omega2, V1)
            dV_TOT = dV_INC + dV_AOP + dV_SMA + dV_RAAN
        else:
            dV_TOT = dV_INC + dV_AOP + dV_SMA

        return dV_TOT

def dV_matrix_generation(debris_data, RAAN_maneuver = False):
    """ 
    Function computing the delta_v matrix of all possible transfer from a debris to another

    Arguments:
        debris_data (DataFrame): Data of the debris considered

        RAAN_maneuver (boolean) - False by default : Boolean indicating if RAAN maneuvers are taken in account in the delta V computation

    Returns:
        dV_matrix (Matrix) : delta-V matrix (triangular superior) recording all possible values of delta_v for the given set of debris
    """


    dV_matrix = np.zeros((len(debris_data), len(debris_data)))

    for i in range(len(debris_data)):
        dV_matrix[i][i]=0
        for j in range(i+1, len(debris_data)):
            dV=total_dV_computation(debris_data, i, j, RAAN_maneuver = RAAN_maneuver)
            dV_matrix[i][j] = dV
            dV_matrix[j][i] = dV

    return dV_matrix
    

def print_dV_matrix(dV_matrix=dV_matrix_generation(DDL.convertTLEtoDF(DDL.recoveringDebrisData()))):
    
    # Setting the font

    rc('font', **{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)
    dV_max=1.5
    ticks=np.arange(1,constants.N_DEBRIS+1,2)
    extent=(0.5,constants.N_DEBRIS+0.5,0.5,constants.N_DEBRIS+0.5)
    cmap0 = LinearSegmentedColormap.from_list('', ['green', 'white'])
    plt.imshow(dV_matrix,vmin=0, vmax=dV_max,interpolation='nearest',extent=extent, origin='lower', cmap=cmap0)
    plt.xlabel('Debris ID')
    plt.ylabel('Debris ID')
    plt.xticks(ticks)
    plt.yticks(ticks)
    colorbar=plt.colorbar()
    colorbar.set_label('$ \Delta V $ [km/s]')

    plt.show()

if __name__ == '__main__':
    
    # Loading and sorting of the debris data
    debris_data = DDL.convertTLEtoDF(DDL.recoveringDebrisData())
    
    # Computation of the delta-V matrix 
    dV_matrix = dV_matrix_generation(debris_data)
    dV_matrix_raan = dV_matrix_generation(debris_data, RAAN_maneuver = True)
    print_dV_matrix()
