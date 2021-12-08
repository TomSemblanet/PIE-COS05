import numpy as np 
import math as mt

def cart2kep(R, V, mu):
	""" Computes the S/C orbital elements from its cartesian coordinates """

	eps = 1e-10

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
		if e > eps:
			w = np.arccos(np.dot(N, E) / n / e)
			if E[2] < 0:
				w = 2*np.pi - w
		else:
			w = 0
	else:
		w = 0

	if e > eps:
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

	return [a, e, i, w, RA, TA]



def kep2cart(coe, mu) : 

	a, e, i, w, RA, TA = coe

	h = np.sqrt(mu * a * (1 - e**2))

	rp = (h**2/mu)*(1/(1+e*np.cos(TA)))*np.array([np.cos(TA), np.sin(TA), 0])
	vp = mu/h*np.array([-np.sin(TA), e+np.cos(TA), 0])


	R1 = np.array([ [np.cos(w), -np.sin(w), 0.],
			        [np.sin(w),  np.cos(w), 0.],
			        [             0.,               0., 1.] ])

	R2 = np.array([ [1.,               0.,               0.],
		            [0.,      np.cos(i),     -np.sin(i)],
		            [0.,      np.sin(i),      np.cos(i)] ])

	R3 = np.array([ [np.cos(RA), -np.sin(RA),     0.],
		            [np.sin(RA),  np.cos(RA),     0.],
		            [             0.,               0.,     1.] ])

	M = R3.dot(R2.dot(R1))

	r = M.dot(rp)
	v = M.dot(vp)
	
	return [-r, -v]



if __name__ == '__main__':

	r = np.array([42000, 0, 0])
	v = np.array([0, 2, 0])
	mu = 398600.4418

	coe = cart2kep(r, v, mu)
	print(kep2cart(coe, mu))






