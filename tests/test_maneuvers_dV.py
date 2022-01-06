import numpy as np
import unittest

from regroupement.dV_computations.maneuvers_dV import SMA_dV, INC_dV, AOP_dV
from utils.constants import mu_EARTH, R_EARTH

class Test_maneuvers_dV(unittest.TestCase):

    def test_INC_dV(self):
        """ Test of the `INC_dV` function """

        i1, i2 = 20 * np.pi/180, 30 * np.pi/180
        V = np.sqrt(mu_EARTH / (R_EARTH + 1000))

        dV_INC_computed = INC_dV(i1, i2, V)
        dV_INC_analytic = 1.28067

        # Equality verification
        self.assertTrue(dV_INC_analytic == round(dV_INC_computed, 5))


    def test_SMA_dV(self):
        """ Test of the `SMA_dV` function """

        a1, a2 = R_EARTH + 1000, R_EARTH + 2000

        dV_SMA_computed = SMA_dV(a1, a2)
        dV_SMA_analytic = 0.49747

        # Equality verification
        self.assertTrue(dV_SMA_analytic == round(dV_SMA_computed, 5))


    def test_AOP_dV(self):
        """ Test of the `AOP_dV` function """

        a, e, i, W = R_EARTH + 1000, 0.1, 10 * np.pi/180, 30 * np.pi/180
        mass = 300
        
        w1, w2 = 10 * np.pi/180, 30 * np.pi/180

        dV_AOP_computed = AOP_dV(w1, w2, W, a, e, i, mass)
        dV_AOP_analytic = 2.56445

        # Equality verification
        self.assertTrue(dV_AOP_analytic == round(dV_AOP_computed, 5))


if __name__ == '__main__':

    unittest.main()

