""" This module computes intermadiary data like velocity, angular momentum, effect of J2 perturbation on RAAN """

import numpy as np
from numpy import linalg as LA
import utilitaires
from utilitaires.coc import kep2cart


def calcul_V(a):
    # param gravitationnel standard terreste mu = G*mTerre
    mu = 398600.4418
    V = np.sqrt(mu/a)
    return V


def angular_momentum(r_v):
    r = r_v[0]
    v = r_v[1]
    return np.cross(r, v)


def derive_omega(a, i, w, M, e, omega, nu, mass):
    alpha = w + M
    mu = 398600.4418
    orbital_param = (a, e, i, w, omega, nu)
    position_velocity = kep2cart(orbital_param, mu)
    h = angular_momentum(position_velocity)
    n = LA.norm(np.cross(np.array[0, 0, 1], h))
    position = position_velocity[0]
    Fgravi = ((mu*mass)/LA.norm(position)**3)*position
    W = np.dot(Fgravi, h)/LA.norm(h)
    derive = (np.sin(alpha)*W)/(n*a*np.sin(i))
    return derive

# dérive relative du RAAN (effet gravitationnel J2) entre deux orbites


def derive_relative(derive1, derive2):
    return derive1 - derive2
