import numpy as np


def delta_V(mu, r1, r2):
    deltaV_1 = np.sqrt((2*mu*r2/(r1*(r1+r2)))) - np.sqrt(mu/r1)
    deltaV_2 = np.sqrt(mu/r2) - np.sqrt((2*mu*r1/(r2*(r1+r2))))
    return deltaV_1, deltaV_2


def delta_m(deltaV, m0):
    isp = 300
    g0 = 9.80665
    return m0*(1-np.exp(-deltaV*1000/(g0*isp)))


# param gravitationnel standard terreste mu = G*mTerre
mu = 398600.4418

# On suppose qu'on part d'une orbite a = 700km)
r1 = np.array([700, 841, 848, 853])
r2 = np.array([841, 848, 853, 863])

delta_m_total = 0
m_vide = 2250  # kg
m0_ergol = 175
m0 = m_vide + m0_ergol


for i in range(r1.shape[0]):
    print(f"orbit transfer number {i}")
    delta_V_init, delta_V_final = delta_V(mu, r1[i], r2[i])
    delta_V_total = delta_V_init + delta_V_final
    print("delta_V_total en km/s")
    print(delta_V_total)
    delta_m_init = delta_m(delta_V_init, m0)
    print("delta_m_init en kg")
    print(delta_m_init)
    m0 -= delta_m_init
    delta_m_final = delta_m(delta_V_final, m0)
    print("delta_m_final en kg")
    print(delta_m_final)
    delta_m_total += delta_m_init + delta_m_final


print("masse totale d'ergols consommée " + str(delta_m_total) + " kg")
