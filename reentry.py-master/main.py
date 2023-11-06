# Import the necessary functions or modules from reentry.py
from reentry import *

# Define your simulation parameters as a dictionary
sim = {
    "max_it": 1000000,
    "delta_t": 0.1,
    "entry_interface": 120e3,   # QARMAN mission
    "fpa": -8.1,
    "velocity": 7.57e3,
    "stop_alt": 0
}

# Define your spacecraft parameters as a dictionary
craft = {
    "name": "test 29/10/2023",
    "ballistic_coef": 146,
    "lift_drag": 0.24
}


# Define the planet object with a "radius" attribute
planet = Earth()

# Call the sim_run function with the defined parameters


# Run the simulation
x, y, vx, vy, ax, ay, t = sim_run(sim, planet, craft)

title = r'%s ($\beta=%.0f\mathrm{kg}/\mathrm{m}^2$) @ %s' % (craft['name'], craft['ballistic_coef'], planet.name)
label = r'$fpa=%.2f, v_\mathrm{EI}=%.0f\mathrm{ms}^{-1}$' % (sim['fpa'], sim['velocity'])

# convert to useful co-oridnates
lat, alt = planet.polar(x, y)
downrange = lat * planet.radius
vn = np.linalg.norm([vx, vy])

# Compute the velocity magnitude
v = np.sqrt(vx ** 2.0 + vy ** 2.0)

# Get the axial load ((ax,ay) projected onto (vx,vy)) --> so also need to consider weight?
aa = np.abs((ax * vx + ay * vy) / v)

# Time and distance to go
tti = np.max(t) - t
dtg = np.max(downrange) - downrange

f1 = plt.figure(figsize=(6, 2))
do_plot(
        'downrange (km)', downrange / 1e3,
        'altitude (km)', alt / 1e3,
        label, title, '%s-traj.png' % craft['name']

    )

f2 = plt.figure(figsize=(4, 3))
do_plot(
        'axial loads (g)', aa / 9.81,
        'altitude (km)', alt / 1e3,
        label, title, '%s-load_alt.png' % craft['name']

    )

f3 = plt.figure(figsize=(4, 3))
do_plot(
        'time since EI (s)', t,
        'axial loads (g)', aa / 9.81,
        label, title, '%s-load_time.png' % craft['name']
    )

f4 = plt.figure(figsize=(4, 3))
do_plot(
        'distance to splashdown (km)', dtg / 1e3,
        'time to parachute deploy (s)', tti / 60.0,
        label, title, '%s-dtg.png' % craft['name']
    )

f5 = plt.figure(figsize=(4, 3))
do_plot(
        'velocity (km/s)', v / 1e3,
        'altitude (km)', alt / 1e3,
        label, title, '%s-vel.png' % craft['name']
    )

plt.close()

