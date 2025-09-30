import matplotlib.pyplot as plt
from rocket import Rocket

def simulate(rocket, duration=5, dt=0.05):
    x_vals = []
    y_vals = []
    time_elapsed = 0

    max_height = 0
    while time_elapsed < duration:
        rocket.update(dt)
        if rocket.y < 0:
            break
        x_vals.append(rocket.x)
        y_vals.append(rocket.y)
        if rocket.y > max_height:
            max_height = rocket.y
        time_elapsed += dt

    final_distance = rocket.x

    return x_vals, y_vals, max_height, final_distance, time_elapsed

    for _ in range(int(duration/dt)):
        rocket.update(dt)
        if rocket.y < 0:
            break
        x_vals.append(rocket.x)
        y_vals.append(rocket.y)

    return x_vals, y_vals

# Example fin shapes with drag coeffs
fin_shapes = {
    'trapezoidal': 0.4,
    'rectangular': 0.5,
    'elliptical': 0.3,
    'delta': 0.6
}

for shape, drag in fin_shapes.items():
    r = Rocket(mass=1.5, fin_shape=shape, drag_coeff=drag, lift_coeff=0.1)
    x, y, max_h, dist, time_flight = simulate(r)

    # plot with fancy label
    plt.plot(x, y, label=f"{shape} (MaxH={round(max_h,1)}m, Dist={round(dist,1)}m)")

    # print to terminal
    print(f"\nFin: {shape.upper()}")
    print(f"→ Max Height: {max_h:.2f} m")
    print(f"→ Distance: {dist:.2f} m")
    print(f"→ Time of Flight: {time_flight:.2f} s")

# Then show it
plt.xlabel('Distance (m)')
plt.ylabel('Height (m)')
plt.title('Rocket Trajectories by Fin Shape')
plt.legend()
plt.grid(True)
plt.show()

