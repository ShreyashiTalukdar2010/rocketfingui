import matplotlib.pyplot as plt
import matplotlib.animation as animation
from rocket import Rocket

# Simulate flight and return positions
def simulate_flight(rocket, dt=0.05, duration=10):
    positions = []
    time_elapsed = 0
    while time_elapsed < duration:
        rocket.update(dt)
        if rocket.y < 0:  # rocket hit ground
            break
        positions.append((rocket.x, rocket.y))
        time_elapsed += dt
    return positions

# Fin shapes with drag coeffs
fin_shapes = {
    'trapezoidal': 0.4,
    'rectangular': 0.5,
    'elliptical': 0.3,
    'delta': 0.6
}

# Simulate each rocket
all_positions = {}
for shape, drag in fin_shapes.items():
    r = Rocket(mass=1.5, fin_shape=shape, drag_coeff=drag, lift_coeff=0.1)
    all_positions[shape] = simulate_flight(r)

# Figure setup
fig, ax = plt.subplots()
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Height (m)")
ax.set_title("Rocket Flight Animation – Different Fin Shapes")

# Axis limits (find global max)
all_x = [p[0] for pos in all_positions.values() for p in pos]
all_y = [p[1] for pos in all_positions.values() for p in pos]
ax.set_xlim(0, max(all_x) * 1.1)
ax.set_ylim(0, max(all_y) * 1.1)

# Create one dot + path line per rocket
colors = {
    'trapezoidal': 'r',
    'rectangular': 'g',
    'elliptical': 'b',
    'delta': 'm'
}
dots = {}
paths = {}
for shape, pos in all_positions.items():
    dot, = ax.plot([], [], 'o', color=colors[shape], markersize=8, label=shape)
    path, = ax.plot([], [], '-', color=colors[shape], linewidth=1.5)
    dots[shape] = dot
    paths[shape] = path

ax.legend()

# Init function
def init():
    for shape in fin_shapes:
        dots[shape].set_data([], [])
        paths[shape].set_data([], [])
    return list(dots.values()) + list(paths.values())

# Update function
def update(frame):
    for shape, pos in all_positions.items():
        if frame < len(pos):
            x, y = pos[frame]
            dots[shape].set_data([x], [y])
            paths[shape].set_data([p[0] for p in pos[:frame+1]],
                                  [p[1] for p in pos[:frame+1]])
    return list(dots.values()) + list(paths.values())

# Max number of frames (longest flight)
max_frames = max(len(pos) for pos in all_positions.values())

ani = animation.FuncAnimation(fig, update, frames=max_frames,
                              init_func=init, blit=True, interval=50, repeat=False)

plt.show()
