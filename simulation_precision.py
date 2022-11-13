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

dt = 36 # časový krok [s] // lze zmenšit pro přesnější simulaci

# inicializace
t = 0
r = r0
v = v0
data = []

while True:
    data.append([t,r[0],r[1], v[0], v[1]])
    abs_r = np.sqrt(r[0]**2 + r[1]**2)
    a = -G*(M/abs_r**3)*r
    v += a*dt
    r += v*dt
    t += dt
    # ukonči cyklus, až po hodně obězích
    if len(data) > 500000:
        break

tt, xx, yy, vvx, vvy = np.hsplit(np.array(data),5)
xx = xx/AU # konvertovat do [AU]
yy = yy/AU # konvertovat do [AU]

hybnost = m * np.sqrt(vvx**2 + vvy**2)
sun_dist = np.sqrt(xx**2 + yy**2)
moment_hybnosti = sun_dist * hybnost

fig, ax = plt.subplots()

ax.plot(tt, moment_hybnosti, "-", color="red", linewidth=1)
ax.set_title("Změna momentu hybnosti v závislosti na čase")
ax.set_ylabel("moment hybnosti [kg·m2·s−1]")
ax.set_xlabel("čas [s]")
ax.grid()

plt.show()
