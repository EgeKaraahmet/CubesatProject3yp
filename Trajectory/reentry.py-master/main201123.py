# Import the necessary functions or modules from reentry.py
from reentry201123 import *
from matplotlib import pyplot as plt


# Define your simulation parameters as a dictionary
sim = {
    "max_it": 10000000,
    "delta_t": 1,
    "entry_interface": 200e3,   # QARMAN mission
    "fpa": 0,
    "velocity": 7788,          
    "stop_alt": 5e3
}

# Define your spacecraft parameters as a dictionary
craft = {
    "name": "Prototype ",
    "ballistic_coef": 112,
    "lift_drag": 0,
    "ballistic_coef_parachute": 2,
}


# Define the planet object with a "radius" attribute
planet = Earth()


# Run the simulation
x, y, vx, vy, ax, ay, t = sim_run141123(sim, planet, craft)

title = r'%s ($\beta=%.0f\mathrm{kg}/\mathrm{m}^2$) @ %s' % (craft['name'], craft['ballistic_coef'], planet.name)
label = r'$fpa=%.2f, v_\mathrm{EI}=%.0f\mathrm{ms}^{-1}$' % (sim['fpa'], sim['velocity'])

# convert to useful co-oridnates
lat, alt = planet.polar(x, y)
downrange = lat * planet.radius
vn = np.linalg.norm([vx, vy])

# Compute the velocity magnitude
v = np.sqrt(vx ** 2.0 + vy ** 2.0)

# Get the axial load ((ax,ay) projected onto (vx,vy)) --> so also need to consider weight?
# aa = np.abs((ax * vx + ay * vy) / v)

aa = np.sqrt((ax **2 + ay **2))


# Time and distance to go
tti = np.max(t) - t
dtg = np.max(downrange) - downrange


def count_max_positive_occurrences(lat):
    oscillation_count = 0
    max_positive_value = float('-inf')

    for value in lat:
        if value > max_positive_value:
            max_positive_value = value
        elif value < 0:
            # A negative value indicates the end of the current oscillation
            if max_positive_value != float('-inf'):
                oscillation_count += 1
            max_positive_value = float('-inf')

    # Check for the last oscillation
    if max_positive_value != float('-inf'):
        oscillation_count += 1

    return oscillation_count

oscillation_count = count_max_positive_occurrences(lat)
print("Number of oscillations reaching maximum positive value:", oscillation_count)

num_points = len(alt)
downrange_flat_earth = np.linspace(0, (oscillation_count-1) * 2*np.pi* planet.radius + np.pi* planet.radius , num_points)
no_of_Earth_circumference = downrange_flat_earth / (2*np.pi*planet.radius)





# Create a single figure and multiple subplots
fig, axs = plt.subplots(2, 3, figsize=(15, 7.5))  # 2 rows and 3 columns of subplots

# Plot the data on each subplot
do_plot(axs[0, 0], 'Downrange (Mm)', downrange_flat_earth/1e6, 'altitude (km)', alt / 1e3, label, title)
# do_plot(axs[0, 0], 'Downrange (km)', downrange / 1e3, 'altitude (km)', alt / 1e3, label, title)
do_plot(axs[1, 0], 'No. of completed Earth Circumference', no_of_Earth_circumference, 'altitude (km)', alt / 1e3, label, title)
# do_plot(axs[1, 0], 'distance to splashdown (km)', dtg / 1e3, 'time to parachute deploy (s)', tti / 60.0, label, title)

do_plot(axs[0, 1], 'axial loads (g)', aa / 9.81, 'altitude (km)', alt / 1e3, label, title)

# Limiting the axies
ax02 = axs[0, 2]
do_plot(ax02, 'time since EI (s)', (t*86400/(86400)), 'axial loads (g)', aa / 9.81, label, title)
# Set the x-axis limits to start from t/86400 = 3
min_x = 194900        # with earth's rotation
# min_x = 183050          # without earth's rotation
# max_x = max(t/86400)  # Calculate the maximum value dynamically
max_x = 195350        # with earth's rotation
# max_x = 2.123*86400     # without earth's rotation
ax02.set_xlim(min_x, max_x)


do_plot(axs[1, 1], 'velocity (km/s)', v / 1e3, 'altitude (km)', alt / 1e3, label, title)

ax12 = axs[1, 2]
do_plot(ax12,'time (day)', t/(86400),'altitude (km)' , alt / 1e3, label, title)
# do_plot(ax12,'time (s)', t,'altitude (km)' , alt / 1e3, label, title)
# min_x = 277500
# max_x = max(t)
# ax12.set_xlim(min_x,max_x)
# min_y = 0
# max_y = 60
# ax12.set_ylim(min_y, max_y)


# Adjust the layout
plt.tight_layout()

# Show the plots
plt.show()