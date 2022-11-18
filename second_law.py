import numpy as np
from matplotlib.widgets import Slider
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
v0 = np.array([0, 2e4]) # pocatecni vektor rychlosti [m/s] - více excentrická dráha pro lepší ilustraci

dt = 3600 # časový krok [s] // lze zmenšit pro přesnější simulaci

# inicializace
t = 0
r = r0
v = v0
S = None
accumulative_area = 0
counter_cross_x = 0
data = []

while True:
    data.append([t,r[0],r[1], S, accumulative_area])
    abs_r = np.sqrt(r[0]**2+r[1]**2)
    a = -G*(M/abs_r**3)*r
    v += a*dt
    r += v*dt
    t += dt

    if len(data) > 1:
        # získá poslední a předposlední záznam lokace
        a_coor = (data[-1][1], data[-1][2])
        b_coor = (data[-2][1], data[-2][2])
        c_coor = (0, 0) # souřadnice Slunce

        # pomocí Pythagorovy věty dopočítá vzdálenosti mezi třemi body
        AC = np.sqrt(a_coor[0]**2 + a_coor[1]**2)
        BC = np.sqrt(b_coor[0]**2 + b_coor[1]**2)
        AB = np.sqrt(abs(a_coor[0]-b_coor[0])**2 + abs(a_coor[1]-b_coor[1])**2)

        # Heronův vzorec na obsah trojúhělníku z SSS
        s = (AC+BC+AB)/2
        S = np.sqrt(s*(s-AC)*(s-AB)*(s-BC))
        accumulative_area += S

    # počítá, kolikrát se změní znaménko na y-ové ose -> planeta oběhla 180deg, půl oběhu
    if np.sign(r[1]) == np.sign(data[-1][2])*-1 or np.sign(r[1]) == 0:
        counter_cross_x += 1

    # okonči cyklus, pokud planeta dvakrát prošla x-ovou osou
    if counter_cross_x == 2:
        break

print(f"Počet oběhu: {counter_cross_x/2}")

tt, xx, yy, ss, ss_accumulative = np.hsplit(np.array(data), 5)
xx = xx/AU # konvertovat do [AU]
yy = yy/AU # konvertovat do [AU]

"""
Zpracování dat pro ověření druhého Keplrova zákonu,
tedy výpočet plochy opsané průvodičem za určitý čas
"""

def get_swept_area(orbit_percentage):
    orbit_time = int(len(tt) * (orbit_percentage/100))
    a_coor = (xx[orbit_time], yy[orbit_time])
    b_coor = (xx[orbit_time+1], yy[orbit_time+1])
    c_coor = (0, 0)
    S = ss[orbit_time]
    return orbit_time, a_coor, b_coor, c_coor, S

def count_rsu():
    min_area = None
    max_area = None
    for percentage in np.arange(0,100,0.2):
        orbit_time, a_coor, b_coor, c_coor, S = get_swept_area(percentage)
        if min_area == None or min_area > S:
            min_area = S
        if max_area == None or max_area < S:
            max_area = S

    diff = abs(max_area - min_area)
    mid_value = max_area - diff/2
    rsu = (diff/2)/mid_value # relativní směrodatná odchylka
    return rsu[0]

rsu = count_rsu()

def update_area(orbit_percentage):
    for patch in ax[0].patches:
        patch.remove()
    for txt in ax[0].texts:
        txt.remove()

    orbit_time, a_coor, b_coor, c_coor, S = get_swept_area(orbit_percentage)
    area_triangle = plt.Polygon(np.array([[xx[orbit_time], yy[orbit_time]],[xx[orbit_time+1], yy[orbit_time+1]],[0,0]], dtype=object), color="red", label=f"plocha opsaná průvodičem planety za {dt} s")
    ax[0].add_patch(area_triangle)

    if a_coor[0] < b_coor[0]:
        ax[0].text(s="A", x=a_coor[0]-0.03, y=a_coor[1])
        ax[0].text(s="B", x=b_coor[0]+0.03, y=b_coor[1])
    elif a_coor[0] >= b_coor[0]:
        ax[0].text(s="A", x=a_coor[0]+0.03, y=a_coor[1])
        ax[0].text(s="B", x=b_coor[0]-0.03, y=b_coor[1])

    ax[0].text(s="C", x=c_coor[0], y=c_coor[1])
    ax[0].text(x=np.min(xx), y=np.min(yy),
            backgroundcolor="#ffffffcf", fontsize="x-small",
            s=f"""
SOUŘADNICE VYKRESLENÉHO TROJÚHELNÍKU
A = ({a_coor[0][0]:.6e}, {a_coor[1][0]:.6e})
B = ({b_coor[0][0]:.6e}, {b_coor[1][0]:.6e})
C = {c_coor[0], c_coor[1]}\n
S = {S[0]:.12e} AU^2
δ = {rsu:.12e}
""")


"""
Vygenerování simulace opsané plochy s posuvníkem a
vyrobení výsledného grafu pro potvrzení zákonu
"""

fig, ax = plt.subplots(nrows=1, ncols=2)

fig.subplots_adjust(bottom=0.2) # uvolní místo dole pro slider
progress_ax= fig.add_axes([0.2, 0.05, 0.6, 0.05])
progress_slider = Slider(
    ax=progress_ax,
    label="procentuální část oběhu",
    valmin=1,
    valmax=99,
    valinit=50,
    valstep=1.0,
    initcolor="none"
)

update_area(40)
progress_slider.on_changed(update_area)

ax[0].plot(xx,yy,"o", color="blue", markersize=1, label="oběžná dráha planety") # oběžná dráha
ax[0].plot(0,0, "o", color="orange", markersize=10, label="Slunce") # Slunce v jednom ohnisku
ax[0].set_ylabel("y [AU]")
ax[0].set_xlabel("x [AU]")
ax[0].set_title("Simulace pohybu planety a obsahu plochy")
ax[0].axis("scaled")
ax[0].grid()
ax[0].legend(loc="upper right", fontsize="x-small")

ax[1].plot(tt, ss_accumulative, "-", color="red", linewidth=1)
ax[1].set_ylabel("kumulativní opsaná plocha [AU²]")
ax[1].set_xlabel("čas [s]")
ax[1].grid()
ax[1].set_title("Kumulativní plocha v závislosti na čase")


plt.suptitle("Potvrzení druhého Keplerova zákona")
plt.show()
