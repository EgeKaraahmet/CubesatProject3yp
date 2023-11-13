# Import the necessary functions or modules from reentry.py
from reentry import *
from matplotlib import pyplot as plt


# Define your simulation parameters as a dictionary
sim = {
    "max_it": 100000,
    "delta_t": 0.1,
    "entry_interface": 500e3,   # QARMAN mission
    "fpa": 0,
    "velocity": 7500,
    "stop_alt": 0
}

# Define your spacecraft parameters as a dictionary
craft = {
    "name": "Prototype ",
    "ballistic_coef": 112,
    "lift_drag": 0
}


# Define the planet object with a "radius" attribute
planet = Earth()


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






# Create a single figure and multiple subplots
fig, axs = plt.subplots(2, 3, figsize=(15, 7.5))  # 2 rows and 3 columns of subplots

# Plot the data on each subplot
do_plot(axs[0, 0], 'downrange (km)', downrange / 1e3, 'altitude (km)', alt / 1e3, label, title)
do_plot(axs[0, 1], 'axial loads (g)', aa / 9.81, 'altitude (km)', alt / 1e3, label, title)
do_plot(axs[0, 2], 'time since EI (s)', t, 'axial loads (g)', aa / 9.81, label, title)
do_plot(axs[1, 0], 'distance to splashdown (km)', dtg / 1e3, 'time to parachute deploy (s)', tti / 60.0, label, title)
do_plot(axs[1, 1], 'velocity (km/s)', v / 1e3, 'altitude (km)', alt / 1e3, label, title)
do_plot(axs[1, 2],'altitude (km)' , alt / 1e3,'latitude (deg)', np.degrees(lat), label, title)

# Adjust the layout
plt.tight_layout()

# Show the plots
plt.show()

