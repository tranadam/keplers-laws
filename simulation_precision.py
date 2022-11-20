import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

"""
Výpočet pohybu planety za vlivu pouze gravitační síly Slunce vycházeje
z hmotností těles, počáteční rychlosti a počáteční vzdálenosti
"""

m = 5.972e24 # hmotnost planety Země [kg]
M = 1.989e30 # hmotnost Slunce [kg]
G = constants.G # gravitační konstanta [m^3 * kg^-1 * s^-2]
AU = constants.astronomical_unit # astronomická jednotka [m]
r0 = np.array([-1.0167*AU, 0*AU]) # počáteční vektor polohy - vzdálenost Země v aféliu
v0 = np.array([0, 2e4]) # počáteční vektor rychlosti [m/s] - více excentrická dráha pro lepší ilustraci

dt = 3600 # časový krok [s] // lze zmenšit pro přesnější simulaci

# inicializace
t = 0
r = r0
v = v0
counter_cross_x = 0
data = []
periheliums = []

while True:
    data.append([t,r[0],r[1], v[0], v[1]])
    abs_r = np.sqrt(r[0]**2 + r[1]**2)
    a = -G*(M/abs_r**3)*r
    v += a*dt
    r += v*dt
    t += dt

    # počítá, kolikrát se změní znaménko na y-ové ose -> planeta oběhla 180deg, půl oběhu
    if np.sign(r[1]) == np.sign(data[-1][2])*-1 or np.sign(r[1]) == 0:
        counter_cross_x += 1

        # pokud je počítadlo liché, náchází se planeta v perihéliu
        if counter_cross_x%2 == 1:
            # spočítá vzdálenost od Slunce, převede na AU a zaokrouhlí
            perihelium_dist = float(f"{(np.sqrt(r[0]**2 + r[1]**2)/AU):.4f}")
            periheliums.append(perihelium_dist)

    # okonči cyklus, pokud planeta prošla x-ovou osou dostkrát
    if counter_cross_x == 200:
        break

tt, xx, yy, vvx, vvy = np.hsplit(np.array(data),5)

num_of_orbits = counter_cross_x/2
orbit_time = t/(counter_cross_x/2) # v sekundách
orbit_days = orbit_time/60/60/24
print(f"Počet oběhu: {num_of_orbits}\nČas oběhu {orbit_days}")

sli = slice(0,-1,int(len(xx)/100)) # vyber pouze 100 záznamů
xx_100 = xx[sli]
yy_100 = yy[sli]
vvx_100 = vvx[sli]
vvy_100 = vvy[sli]
tt_100 = tt[sli]

rr_100 = np.sqrt(xx_100**2 + yy_100**2)
vv_100 = np.sqrt(vvx_100**2+vvy_100**2)

e_potential = -(G*M*m/rr_100)
e_mechanical = 1/2*m*vv_100**2
total_energy = (e_mechanical + e_potential) / 10**6 # to MJ

fig, ax = plt.subplots(nrows=1, ncols=2)

ax[0].plot(tt_100/orbit_time, total_energy, "-", color="red", linewidth=1)
ax[0].set_title("Vývoj celkové energie v čase")
ax[0].set_ylabel("energie [MJ]")
ax[0].set_xlabel("počet oběhů")
ax[0].grid()

ax[1].plot(np.arange(0,len(periheliums)), periheliums, "-", color="red",linewidth=1)
ax[1].set_title("Změna vzdálenosti planety od Slunce v perihéliu")
ax[1].set_ylabel("vzdálenost [AU]")
ax[1].set_xlabel("počet oběhů")
ax[1].grid()

plt.show()
