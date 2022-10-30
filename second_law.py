import numpy as np
from matplotlib.widgets import Slider
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

m = 5.972e24 # hmotnost planety Země [kg]
M = 1.989e30 # hmotnost Slunce [kg]
G = 6.6743e-11 # gravitační konstanta [m^3 * kg^-1 * s^-2]
AU = 149597870700.0 # astronomická jednotka [m]
r0 = np.array([-1.0167*AU, 0*AU]) # počáteční vektor polohy - vzdálenost Země v aféliu
v0 = np.array([0, 2e4]) # pocatecni vektor rychlosti [m/s] - více excentrická dráha
# v0 = np.array([0, 29290.0]) # pocatecni vektor rychlosti [m/s] - rychlost Země v aféliu

dt = 360 # časový krok [s] // lze zmenšit pro přesnější simulaci

# inicializace
t = 0
r = r0
v = v0
data = []

while True:
    abs_r = np.sqrt(r[0]**2+r[1]**2)
    a = -G*(M/abs_r**3)*r
    data.append([t,r[0],r[1],v[0],v[1]])
    v = v + a*dt
    r = r + v*dt
    t += dt
    # okonči cyklus, pokud byla oběhnuta celá dráha jednou // momentální vektorová poloha je mezi úplně prvním a druhým záznamem
    if (t > dt*2) and (r[0] >= data[0][1] and r[0] < data[1][1]) and (r[1] >= data[0][2] and r[1] < data[1][2]):
        break
tt, xx, yy, vvx, vvy = np.hsplit(np.array(data),5)
xx = xx/AU # [AU]
yy = yy/AU # [AU]
vvx = vvx/AU # [AU/s]
vvy = vvy/AU # [AU/s]

ax.plot(xx,yy,"o", color="blue", markersize=1, label="oběžná dráha planety") # oběžná dráha
ax.plot(0,0, "o", color="orange", markersize=10, label="Slunce") # Slunce v jednom ohnisku

fig.subplots_adjust(bottom=0.2) # uvolní místo dole pro slider
progress_ax= fig.add_axes([0.2, 0.05, 0.6, 0.05])
progress_slider = Slider(
    ax=progress_ax,
    label="Procentuální část oběhu",
    valmin=0,
    valmax=99,
    valinit=50,
    valstep=1.0,
    initcolor="none"
)

def count_swept_area(orbit_percentage):
    orbit_time = int(len(tt) * (orbit_percentage/100))
    a_coor = np.array([xx[orbit_time], yy[orbit_time]])
    b_coor = np.array([xx[orbit_time+1], yy[orbit_time+1]])
    c_coor = np.array([[0],[0]])

    AC = np.sqrt(a_coor[0]**2+a_coor[1]**2)
    BC = np.sqrt(b_coor[0]**2+b_coor[1]**2)
    AB = np.sqrt(abs(a_coor[0]-b_coor[0])**2+abs(a_coor[1]-b_coor[1])**2)

    # Heronův vzorec na obsah trojúhělníku z SSS
    s = (AC+BC+AB)/2
    S = np.sqrt(s*(s-AC)*(s-AB)*(s-BC))
    return orbit_time, a_coor, b_coor, c_coor, S

def count_rsu():
    min_area = None
    max_area = None
    for percentage in np.arange(0,100,0.2):
        orbit_time, a_coor, b_coor, c_coor, S = count_swept_area(percentage)
        if min_area == None or min_area > S:
            min_area = S
        if max_area == None or max_area < S:
            max_area = S

    diff = abs(max_area-min_area)
    mid_value = max_area - diff/2
    rsu = (diff/2)/mid_value # relative standard uncertainty
    return rsu[0]

def update_area(orbit_percentage):
    orbit_time, a_coor, b_coor, c_coor, S = count_swept_area(orbit_percentage)
    for patch in ax.patches:
        patch.remove()
    for txt in ax.texts:
        txt.remove()

    area_triangle = plt.Polygon(np.array([[xx[orbit_time], yy[orbit_time]],[xx[orbit_time+1], yy[orbit_time+1]],[0,0]], dtype=object), color="red", label=f"plocha opsaná průvodičem planety za {dt}s")
    ax.add_patch(area_triangle)
    if a_coor[0] < b_coor[0]:
        ax.text(s="A", x=a_coor[0]-0.03, y=a_coor[1])
        ax.text(s="B", x=b_coor[0]+0.03, y=b_coor[1])
    elif a_coor[0] >= b_coor[0]:
        ax.text(s="A", x=a_coor[0]+0.03, y=a_coor[1])
        ax.text(s="B", x=b_coor[0]-0.03, y=b_coor[1])
    ax.text(s="C", x=c_coor[0], y=c_coor[1])
    ax.text(x=np.amin(xx), y=np.amin(yy),
            backgroundcolor="#ffffffcf", fontsize="x-small",
            s=f"""
SOUŘADNICE VYKRESLENÉHO TROJÚHELNÍKU
A = {a_coor[0][0], a_coor[1][0]}
B = {b_coor[0][0], b_coor[1][0]}
C = {c_coor[0][0], c_coor[1][0]}\n
S = {S[0]} m^2
δ = {count_rsu()}
""")

update_area(50)
progress_slider.on_changed(update_area)

ax.set_ylabel("y [AU]")
ax.set_xlabel("x [AU]")
ax.set_title("Potvrzení druhého Keplerova zákona")
ax.axis("scaled")
ax.grid()
ax.legend(loc="upper right", fontsize="x-small")
plt.show()
