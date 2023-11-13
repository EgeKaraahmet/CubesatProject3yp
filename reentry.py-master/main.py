# Import the necessary functions or modules from reentry.py
from reentry import *
from matplotlib import pyplot as plt


# Define your simulation parameters as a dictionary
sim = {
    "max_it": 10000000,
    "delta_t": 1,
    "entry_interface": 200e3,   # QARMAN mission
    "fpa": 0,
    "velocity": 7788,
    "stop_alt": 4e3
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



num_points = len(alt)
downrange_flat_earth = np.linspace(0, num_points, num_points)
no_of_Earth_circumference = downrange_flat_earth /(2*np.pi*planet.radius)


# Create a single figure and multiple subplots
fig, axs = plt.subplots(2, 3, figsize=(15, 7.5))  # 2 rows and 3 columns of subplots

# Plot the data on each subplot
do_plot(axs[0, 0], 'Downrange (km)', downrange_flat_earth / 1e3, 'altitude (km)', alt / 1e3, label, title)
# do_plot(axs[0, 0], 'Downrange (km)', downrange / 1e3, 'altitude (km)', alt / 1e3, label, title)
do_plot(axs[1, 0], 'No. of completed Earth Circumference', no_of_Earth_circumference, 'altitude (km)', alt / 1e3, label, title)
# do_plot(axs[1, 0], 'distance to splashdown (km)', dtg / 1e3, 'time to parachute deploy (s)', tti / 60.0, label, title)

do_plot(axs[0, 1], 'axial loads (g)', aa / 9.81, 'altitude (km)', alt / 1e3, label, title)

# Limiting the axies
ax = axs[0, 2]
do_plot(ax, 'time since EI (days)', (t/(86400)), 'axial loads (g)', aa / 9.81, label, title)
# Set the x-axis limits to start from t/86400 = 3
min_x = 3.2
max_x = max(t/86400)  # Calculate the maximum value dynamically
ax.set_xlim(min_x, max_x)


do_plot(axs[1, 1], 'velocity (km/s)', v / 1e3, 'altitude (km)', alt / 1e3, label, title)
do_plot(axs[1, 2],'time (days)', t/(86400),'altitude (km)' , alt / 1e3, label, title)

# Adjust the layout
plt.tight_layout()

# Show the plots
plt.show()

