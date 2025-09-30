import math
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from rocket import Rocket  # your rocket class


# --------- Physics helpers ---------
def simulate_flight(rocket, angle_deg=60, speed=50, dt=0.05, duration=10.0):
    """Return list of (x, y) positions until ground or duration reached."""
    positions = []
    t = 0.0

    # set initial velocity here
    angle = math.radians(angle_deg)
    rocket.vx = speed * math.cos(angle)
    rocket.vy = speed * math.sin(angle)

    while t < duration:
        rocket.update(dt)
        if rocket.y < 0:
            break
        positions.append((rocket.x, rocket.y))
        t += dt
    return positions, t

def flight_stats(positions, dt):
    if not positions:
        return 0.0, 0.0, 0.0
    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    max_h = max(ys)
    dist = xs[-1]
    time = len(positions) * dt
    return max_h, dist, time


# --------- GUI ---------
class RocketGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Rocket Fin Simulator 🚀")

        # --- Controls frame ---
        self.controls = ttk.Frame(root, padding=10)
        self.controls.grid(row=0, column=0, sticky="nsw")

        # --- Plot frame ---
        self.plot_frame = ttk.Frame(root, padding=10)
        self.plot_frame.grid(row=0, column=1, sticky="nsew")
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)

        # --- Sliders / Variables ---
        self.fin_var = tk.StringVar(value="elliptical")
        self.mass_var = tk.DoubleVar(value=1.5)
        self.drag_var = tk.DoubleVar(value=0.3)
        self.lift_var = tk.DoubleVar(value=0.10)
        self.angle_var = tk.DoubleVar(value=75.0)
        self.speed_var = tk.DoubleVar(value=50.0)
        self.dt_var = tk.DoubleVar(value=0.05)

        # Fin selector
        ttk.Label(self.controls, text="Fin shape").grid(row=0, column=0, sticky="w")
        self.fin_combo = ttk.Combobox(
            self.controls,
            textvariable=self.fin_var,
            values=["trapezoidal", "rectangular", "elliptical", "delta"],
            state="readonly",
            width=14
        )
        self.fin_combo.grid(row=0, column=1, padx=6, pady=4)

        # Mass / drag / lift / angle / speed / dt sliders
        self._slider("Mass (kg)", self.mass_var, 0.5, 5.0, 0.1, 1)
        self._slider("Drag coeff", self.drag_var, 0.1, 1.0, 0.01, 2)
        self._slider("Lift coeff", self.lift_var, 0.0, 0.5, 0.01, 3)
        self._slider("Angle (deg)", self.angle_var, 15, 90, 1, 4)
        self._slider("Speed (m/s)", self.speed_var, 10, 120, 1, 5)
        self._slider("dt (s)", self.dt_var, 0.01, 0.1, 0.005, 6)

        # --- Buttons ---
        btns = ttk.Frame(self.controls)
        btns.grid(row=7, column=0, columnspan=2, pady=(10, 6), sticky="we")
        self.launch_btn = ttk.Button(btns, text="Launch One 🚀", command=self.launch_one)
        self.launch_btn.grid(row=0, column=0, padx=4)
        self.race_btn = ttk.Button(btns, text="Race All Fins 🏁", command=self.race_all)
        self.race_btn.grid(row=0, column=1, padx=4)
        self.clear_btn = ttk.Button(btns, text="Clear", command=self.clear_plot)
        self.clear_btn.grid(row=0, column=2, padx=4)

        # Mission log
        ttk.Label(self.controls, text="Mission Log").grid(row=8, column=0, columnspan=2, sticky="w", pady=(8, 2))
        self.log = tk.Text(self.controls, width=36, height=10)
        self.log.grid(row=9, column=0, columnspan=2, sticky="we")

        # --- Matplotlib figure ---
        self.fig, self.ax = plt.subplots(figsize=(6.5, 4.5), dpi=100)
        self.ax.set_xlabel("Distance (m)")
        self.ax.set_ylabel("Height (m)")
        self.ax.set_title("Launch view")
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # --- Animations / drawables ---
        self.active_lines = []
        self.active_dots = []
        self.animations = []

        # --- Colors ---
        self.colors = {
            "trapezoidal": "tab:red",
            "rectangular": "tab:green",
            "elliptical": "tab:blue",
            "delta": "tab:purple"
        }

    def _slider(self, label, var, frm, to, res, row):
        ttk.Label(self.controls, text=label).grid(row=row, column=0, sticky="w")
        s = ttk.Scale(self.controls, from_=frm, to=to, orient="horizontal",
                      variable=var)
        s.grid(row=row, column=1, sticky="we", padx=6, pady=2)
        val = ttk.Label(self.controls, textvariable=var)
        val.grid(row=row, column=2, sticky="w")

    # --------- Actions ---------
    def clear_plot(self):
        for line in self.active_lines:
            line.remove()
        for dot in self.active_dots:
            dot.remove()
        self.active_lines.clear()
        self.active_dots.clear()
        self.animations.clear()
        self.ax.cla()
        self.ax.set_xlabel("Distance (m)")
        self.ax.set_ylabel("Height (m)")
        self.ax.set_title("Launch view")
        self.ax.grid(True)
        self.canvas.draw_idle()

    def build_rocket(self, fin_shape, mass, drag, lift, angle_deg, speed_mag):
        r = Rocket(mass=mass, fin_shape=fin_shape, drag_coeff=drag, lift_coeff=lift)
        ang = math.radians(angle_deg)
        r.vx = speed_mag * math.cos(ang)
        r.vy = speed_mag * math.sin(ang)
        return r

    def animate_positions(self, positions, color, label, dt):
        if not positions:
            messagebox.showwarning("No flight", "Rocket didn’t get off the ground!")
            return

        xs = [p[0] for p in positions]
        ys = [p[1] for p in positions]

        self.ax.set_xlim(0, max(xs)*1.1)
        self.ax.set_ylim(0, max(ys)*1.1)

        (line,) = self.ax.plot([], [], '-', linewidth=2, color=color, label=label)
        (dot,) = self.ax.plot([], [], 'o', color=color, markersize=6)
        self.active_lines.append(line)
        self.active_dots.append(dot)
        self.ax.legend(loc="upper right")

        def init():
            line.set_data([], [])
            dot.set_data([], [])
            return line, dot

        def update(frame):
            line.set_data(xs[:frame+1], ys[:frame+1])
            dot.set_data([xs[frame]], [ys[frame]])
            return line, dot

        ani = animation.FuncAnimation(
            self.fig, update, frames=len(xs),
            init_func=init, blit=True, interval=max(10, int(dt*1000)), repeat=False
        )
        self.animations.append(ani)
        self.canvas.draw_idle()

        max_h, dist, time = flight_stats(positions, dt)
        self.log.insert("end", f"{label}: MaxH={max_h:.1f} m | Dist={dist:.1f} m | Time={time:.2f} s\n")
        self.log.see("end")

    def launch_one(self):
        fin = self.fin_var.get()
        mass = self.mass_var.get()
        drag = self.drag_var.get()
        lift = self.lift_var.get()
        angle = self.angle_var.get()
        speed = self.speed_var.get()
        dt = self.dt_var.get()

        rocket = self.build_rocket(fin, mass, drag, lift, angle, speed)
        positions, _t = simulate_flight(rocket, dt=dt, duration=12.0)
        self.animate_positions(positions, self.colors.get(fin, "black"), fin, dt)
    def race_all(self):
        self.clear_plot()
        for shape, drag in fin_shapes.items():
            rocket = Rocket(fin_drag=drag)  # no angle here
            traj, _ = simulate_flight(rocket, angle_deg=60, speed=50)  # pass angle+speed here
            x = [p[0] for p in traj]
            y = [p[1] for p in traj]
            self.ax.plot(x, y, label=shape)
        self.ax.legend()
        self.canvas.draw()
        
# --------- Main ---------
def main():
    root = tk.Tk()
    app = RocketGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()

