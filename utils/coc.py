import numpy as np

from utils.constants import mu_EARTH


def kep2cart(coe, mu):
    """ Converts coordinates of a body from orbital elements (coe) into cartesian coordinates 
    Arguments : 
    
        coe (array) : Orbital elements of the body (SMA, ECC, INC, AOP, RAAN, TA)

        mu (float) : Characteristic parameter of the central body 

    Returns : 
        r (array) : Concatenation of the body's position and velocity (X, Y, Z, VX, VY, VZ) in the geocentric frame
    """

    a, e, i, w, RA, TA = coe

    h = np.sqrt(mu * a * (1 - e**2))

    rp = (h**2/mu)*(1/(1+e*np.cos(TA)))*np.array([np.cos(TA), np.sin(TA), 0])
    vp = mu/h*np.array([-np.sin(TA), e+np.cos(TA), 0])

    R1 = np.array([[np.cos(w), -np.sin(w), 0.],
                   [np.sin(w),  np.cos(w), 0.],
                   [0.,               0., 1.]])

    R2 = np.array([[1.,               0.,               0.],
                   [0.,      np.cos(i),     -np.sin(i)],
                   [0.,      np.sin(i),      np.cos(i)]])

    R3 = np.array([[np.cos(RA), -np.sin(RA),     0.],
                   [np.sin(RA),  np.cos(RA),     0.],
                   [0.,               0.,     1.]])

    M = R3.dot(R2.dot(R1))

    r = M.dot(rp)
    v = M.dot(vp)

    return [r, v]


def cart2kep(R, V, mu):
    """ Converts coordinates of a body from its cartesian coordinates into its orbitals elements (coe)
    
    Arguments : 
    
        R (array) : Position of the body in the ECI frame

        V (array) : Velocity of the body in the ECI frame

        mu (float) : Characteristic parameter of the central body 

    Returns : 
        coe (array) : Body's orbital elements (SMA, ECC, INC, AOP, RAAN, TA)
    """

    r = np.linalg.norm(R)
    v = np.linalg.norm(V)

    vr = np.dot(R, V) / r

    H = np.cross(R, V)
    h = np.linalg.norm(H)

    i = np.arccos(H[2] / h)

    N = np.cross([0, 0, 1], H)
    n = np.linalg.norm(N)

    if n != 0:
        RA = np.arccos(N[0]/n)
        if N[1] < 0:
            RA = 2*np.pi - RA
    else:
        RA = 0

    E = 1/mu * ((v**2 - mu/r)*R - r*vr*V)
    e = np.linalg.norm(E)

    if n != 0:
        if e > 1e-10:
            w = np.arccos(np.dot(N, E) / n / e)
            if E[2] < 0:
                w = 2*np.pi - w
        else:
            w = 0
    else:
        w = 0

    if e > 1e-10:
        TA = np.arccos(np.dot(E, R) / e / r)
        if vr < 0:
            TA = 2*np.pi - TA
    else:
        cp = np.cross(N, R)
        if cp[2] >= 0:
            TA = np.arccos(np.dot(N, R) / n / r)
        else:
            TA = 2*np.pi - np.arccos(np.dot(N, R) / n / r)

    a = h**2 / mu / (1 - e**2)

    coe = [a, e, i, w, RA, TA]

    return coe
