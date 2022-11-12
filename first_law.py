import numpy as np
from matplotlib.widgets import Slider
import matplotlib.pyplot as plt

"""
Výpočet pohybu planety za vlivu pouze gravitační síly Slunce vycházeje
z hmotností těles, počáteční rychlosti a počáteční vzdálenosti
"""

m = 5.972e24 # hmotnost planety Země [kg]
M = 1.989e30 # hmotnost Slunce [kg]
G = 6.6743e-11 # gravitační konstanta [m^3 * kg^-1 * s^-2]
AU = 149597870700.0 # astronomická jednotka [m]
r0 = np.array([-1.0167*AU, 0*AU]) # počáteční vektor polohy - vzdálenost Země v aféliu
# v0 = np.array([0, 29290.0]) # pocatecni vektor rychlosti [m/s] - rychlost Země v aféliu
v0 = np.array([0, 2e4]) # pocatecni vektor rychlosti [m/s] - více excentrická dráha pro lepší ilustraci

dt = 3600 # časový krok [s] // lze zmenšit pro přesnější simulaci

# inicializace
t = 0
r = r0
v = v0
data = []

while True:
    data.append([t,r[0],r[1]])
    abs_r = np.sqrt(r[0]**2 + r[1]**2)
    a = -G*(M/abs_r**3)*r
    v += a*dt
    r += v*dt
    t += dt
    # okonči cyklus, pokud byla oběhnuta celá dráha jednou // momentální vektorová poloha je mezi úplně prvním a druhým záznamem
    if (t > dt*2) and (r[0] >= data[0][1] and r[0] < data[1][1]) and (r[1] >= data[0][2] and r[1] < data[1][2]):
        break
    # ukonči cyklus, pokud simulace trvá moc dlouho // oběh se neuzavírá, nebo je tam až příliš kroků
    if len(data) > 500000:
        break

"""
Zpracování dat pro ověření prvního Keplrova zákonu,
tedy výpočet součtů vzdáleností od ohnisek.
"""

tt, xx, yy = np.hsplit(np.array(data),3)
xx = xx/AU # konvertovat do [AU]
yy = yy/AU # konvertovat do [AU]

def find_ellipse_extremes():
    x_min = None
    x_max = None
    y_min = None
    y_max = None
    for x in xx:
        if x_min == None or x_min > x:
            x_min = x
        if x_max == None or x_max < x:
            x_max = x
    for y in yy:
        if y_min == None or y_min > y:
            y_min = y
        if y_max == None or y_max < y:
            y_max = y
    return x_min, x_max, y_min, y_max

def find_second_focus(extremes_data):
    x_min, x_max, y_min, y_max = extremes_data
    x_focus = x_min + x_max
    y_focus = y_min + y_max
    return x_focus, y_focus

def count_foci_point_distance(orbit_time):
    x_focus = second_focus_coor[0]
    y_focus = second_focus_coor[1]
    x_sun = np.array([0])
    y_sun = np.array([0])
    x_point = xx[orbit_time]
    y_point = yy[orbit_time]

    sun_point_dist = np.sqrt(x_point**2 + y_point**2)
    focus_point_dist = np.sqrt(abs(x_point-x_focus)**2 + abs(y_point-y_focus)**2)
    whole_dist = sun_point_dist + focus_point_dist
    return whole_dist, x_focus, y_focus, x_sun, y_sun, x_point, y_point

def count_rsu():
    min_dist = None
    max_dist = None
    for percentage in np.arange(0,100,0.2):
        orbit_time = int(len(tt) * (percentage/100))
        whole_dist, x_focus, y_focus, x_sun, y_sun, x_point, y_point = count_foci_point_distance(orbit_time)
        if min_dist == None or min_dist > whole_dist:
            min_dist = whole_dist
        if max_dist == None or max_dist < whole_dist:
            max_dist = whole_dist

    diff = abs(max_dist - min_dist)
    mid_value = max_dist - diff/2
    rsu = (diff/2)/mid_value # relativní směrodatná odchylka
    return rsu[0]

def update_plot(orbit_percentage):
    for patch in ax[0].patches:
        patch.remove()
    for txt in ax[0].texts:
        txt.remove()

    orbit_time = int(len(tt) * (orbit_percentage/100))
    whole_dist, x_focus, y_focus, x_sun, y_sun, x_point, y_point = count_foci_point_distance(orbit_time)
    connection_line = plt.Polygon(np.array([[x_focus[0], y_focus[0]],[x_point[0], y_point[0]], [x_sun, y_sun]], dtype=object), color="red",label="spojnice ohniska s bodem na elipse a druhým ohniskem", closed=False, fill=False)
    ax[0].add_patch(connection_line)

    ax[0].text(x=np.amin(xx), y=np.amin(yy),
            backgroundcolor="#ffffffcf", fontsize="x-small",
            s=f"""
SOUŘADNICE BODŮ
A = ({x_point[0]:.6e}, {y_point[0]:.6e})
F₁ = ({x_focus[0]:.6e}, {y_focus[0]:.6e})
F₂, Slunce = {x_sun[0], y_sun[0]}

|AF₁| + |AF₂| = {whole_dist[0]:.8f} AU
δ = {rsu:.8e}
""")

second_focus_coor = (find_second_focus(find_ellipse_extremes()))
rsu = count_rsu()

"""
Vygenerování simulace součtu vzdáleností od ohnisek s posuvníkem
a vyrobení výsledného grafu pro potvrzení zákonu
"""

fig, ax = plt.subplots(nrows=1, ncols=2)

fig.subplots_adjust(bottom=0.2) # uvolní místo dole pro slider
progress_ax= fig.add_axes([0.2, 0.05, 0.6, 0.05])
progress_slider = Slider(
    ax=progress_ax,
    label="Procentuální část oběhu",
    valmin=1,
    valmax=99,
    valinit=40,
    valstep=1.0,
    initcolor="none"
)

update_plot(40)
progress_slider.on_changed(update_plot)

graph_data = []
for index, time in enumerate(tt):
    dist = count_foci_point_distance(index)[0]
    graph_data.append([float(f"{dist[0]:.4f}"), time[0]])
graph_dist, graph_time = np.hsplit(np.array(graph_data), 2)


fig.suptitle("Potvrzení prvního Keplerova zákona")

x_focus = find_second_focus(find_ellipse_extremes())[0]
y_focus = find_second_focus(find_ellipse_extremes())[1]
ax[0].plot(xx,yy,"o", color="blue", markersize=1, label="oběžná dráha planety")
ax[0].plot(0,0, "o", color="orange", markersize=10, label="Slunce")
ax[0].plot(x_focus, y_focus, "o", color="black", markersize=10, label="Druhé ohnisko")
ax[0].set_ylabel("y [AU]")
ax[0].set_xlabel("x [AU]")
ax[0].set_title("Animace pohybu planety a součtu vzdáleností od ohnisek")
ax[0].axis("scaled")
ax[0].grid()
ax[0].legend(loc="upper right", fontsize="x-small")

ax[1].plot(graph_time,graph_dist, "-", color="red", linewidth=1)
ax[1].set_ylabel("vzdálenost od ohnisek [AU]") # čtyři desetinná místa
ax[1].set_xlabel("čas [s]")
ax[1].set_title("vzdálenost od ohnisek v závislosti na čas")
ax[1].grid()


plt.show()
