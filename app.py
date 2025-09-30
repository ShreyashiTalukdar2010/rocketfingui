# app.py — Rocket Fin Simulator (MIT-ready, single-file)
# Requirements: Python 3.10+, matplotlib, tkinter (comes with Python on Windows)
# Run: python app.py

import math
import random
import tkinter as tk
from tkinter import ttk, messagebox

import matplotlib
matplotlib.use("TkAgg")  # Tkinter-friendly backend
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# -----------------------------
# Physics Model (lightweight)
# -----------------------------
class Rocket:
    """
    Simple 2D point-mass rocket with drag; 'fin shape' affects Cd (drag coefficient).
    No thrust stage here — it's a toss/launch demo to study fin-shape effects on drag.
    """
    g = 9.81  # m/s^2

    def __init__(self, mass: float, cd: float, vx0: float, vy0: float):
        self.mass = mass
        self.cd = cd            # drag coefficient (shape-dependent simplification)
        self.x = 0.0
        self.y = 0.0
        self.vx = vx0
        self.vy = vy0

    def step(self, dt: float, wind_ax: float = 0.0):
        """
        Integrate one time step with quadratic drag and gravity.
        wind_ax: horizontal acceleration from wind (+ right, - left).
        """
        # Speed and unit vector
        v = math.hypot(self.vx, self.vy)
        if v > 1e-9:
            ux, uy = self.vx / v, self.vy / v
        else:
            ux, uy = 0.0, 0.0

        # Drag magnitude ~ 0.5 * Cd * v^2 (we absorb rho*A into Cd for simplicity)
        Fd = 0.5 * self.cd * v * v
        Fx_drag = -Fd * ux
        Fy_drag = -Fd * uy

        # Gravity
        Fy_grav = -self.mass * self.g

        # Total forces (add simple wind as horizontal accel * m)
        Fx = Fx_drag + self.mass * wind_ax
        Fy = Fy_drag + Fy_grav

        # Accelerations
        ax = Fx / self.mass
        ay = Fy / self.mass

        # Integrate velocity
        self.vx += ax * dt
        self.vy += ay * dt

        # Integrate position
        self.x += self.vx * dt
        self.y += self.vy * dt


def simulate_flight(mass: float,
                    cd: float,
                    angle_deg: float,
                    speed: float,
                    dt: float = 0.02,
                    duration: float = 12.0,
                    wind_on: bool = False):
    """
    Simulate one rocket and return (positions, stats).
    positions: list[(x,y)] until y < 0 or duration reached.
    stats: dict with max_height, distance, time_of_flight.
    """
    # Initial velocity from angle & speed
    ang = math.radians(angle_deg)
    vx0 = speed * math.cos(ang)
    vy0 = speed * math.sin(ang)

    r = Rocket(mass=mass, cd=cd, vx0=vx0, vy0=vy0)
    positions = []
    t = 0.0
    max_h = 0.0

    # Wind model: small random gusts, changes slowly
    gust = 0.0
    gust_timer = 0.0

    while t < duration:
        if wind_on:
            gust_timer -= dt
            if gust_timer <= 0.0:
                # new gust every 0.4–0.8s, accel between -1.2..1.2 m/s^2
                gust = random.uniform(-1.2, 1.2)
                gust_timer = random.uniform(0.4, 0.8)
            wind_ax = gust
        else:
            wind_ax = 0.0

        r.step(dt, wind_ax=wind_ax)

        if r.y < 0.0:
            break

        positions.append((r.x, r.y))
        if r.y > max_h:
            max_h = r.y
        t += dt

    distance = positions[-1][0] if positions else 0.0
    stats = dict(max_height=max_h, distance=distance, time_of_flight=t)
    return positions, stats


# Fin shapes → drag coefficient (relative)
FINS = {
    "rectangular": 0.55,
    "trapezoidal": 0.48,
    "elliptical": 0.38,
    "delta": 0.62,
}

FIN_COLORS = {
    "rectangular": "tab:blue",
    "trapezoidal": "tab:orange",
    "elliptical": "tab:green",
    "delta": "tab:red",
}


# -----------------------------
# GUI Application
# -----------------------------
class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Rocket Fin Simulator — MIT Maker Portfolio")

        # State
        self.history = []       # list of dicts: {'name', 'positions', 'stats', 'color'}
        self.last_race = []     # same structure but for last race run
        self.animations = []    # keep references to FuncAnimation
        self.current_mode = tk.StringVar(value="Trajectory")

        # Layout: left controls, right plot
        self._build_layout()
        self._build_controls()
        self._build_plot()

    # ----- Layout helpers -----
    def _build_layout(self):
        self.left = ttk.Frame(self.root, padding=10)
        self.left.grid(row=0, column=0, sticky="nsw")

        self.right = ttk.Frame(self.root, padding=10)
        self.right.grid(row=0, column=1, sticky="nsew")

        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)

    def _build_controls(self):
        # Inputs
        r = 0
        ttk.Label(self.left, text="Fin shape").grid(row=r, column=0, sticky="w")
        self.fin_var = tk.StringVar(value="elliptical")
        ttk.Combobox(
            self.left, textvariable=self.fin_var,
            values=list(FINS.keys()), state="readonly", width=15
        ).grid(row=r, column=1, sticky="w", padx=6, pady=4)

        r += 1
        ttk.Label(self.left, text="Angle (deg)").grid(row=r, column=0, sticky="w")
        self.angle_var = tk.DoubleVar(value=60.0)
        self._slider(row=r, var=self.angle_var, frm=10, to=85, res=1)

        r += 1
        ttk.Label(self.left, text="Speed (m/s)").grid(row=r, column=0, sticky="w")
        self.speed_var = tk.DoubleVar(value=50.0)
        self._slider(row=r, var=self.speed_var, frm=10, to=140, res=1)

        r += 1
        ttk.Label(self.left, text="Mass (kg)").grid(row=r, column=0, sticky="w")
        self.mass_var = tk.DoubleVar(value=1.5)
        self._slider(row=r, var=self.mass_var, frm=0.5, to=5.0, res=0.1)

        r += 1
        self.wind_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(self.left, text="Wind (gusts)", variable=self.wind_var).grid(row=r, column=0, columnspan=2, sticky="w", pady=(6, 10))

        # Buttons
        r += 1
        btns = ttk.Frame(self.left)
        btns.grid(row=r, column=0, columnspan=2, pady=(4, 8), sticky="we")
        ttk.Button(btns, text="Single Launch 🚀", command=self.launch_one).grid(row=0, column=0, padx=4)
        ttk.Button(btns, text="Race Mode 🏁", command=self.race_mode).grid(row=0, column=1, padx=4)
        ttk.Button(btns, text="Replay Ghosts 👻", command=self.replay_ghosts).grid(row=0, column=2, padx=4)

        # Graph mode switcher
        r += 1
        ttk.Label(self.left, text="Graph Mode").grid(row=r, column=0, sticky="w")
        gm = ttk.Combobox(self.left, textvariable=self.current_mode, values=["Trajectory", "Bar Chart"], state="readonly", width=15)
        gm.grid(row=r, column=1, sticky="w", padx=6, pady=4)
        gm.bind("<<ComboboxSelected>>", lambda _evt: self.redraw())

        # Live stats
        r += 1
        ttk.Label(self.left, text="Live Stats").grid(row=r, column=0, columnspan=2, sticky="w", pady=(8, 2))
        self.stats_text = tk.Text(self.left, width=34, height=10)
        self.stats_text.grid(row=r+1, column=0, columnspan=2, sticky="we")

        # Clear
        r += 2
        ttk.Button(self.left, text="Clear Plot", command=self.clear_plot).grid(row=r, column=0, pady=8)
        ttk.Button(self.left, text="Clear History", command=self.clear_history).grid(row=r, column=1, pady=8)

    def _slider(self, row, var, frm, to, res):
        s = ttk.Scale(self.left, from_=frm, to=to, orient="horizontal", variable=var)
        s.grid(row=row, column=1, sticky="we", padx=6)
        # live label
        val = ttk.Label(self.left, textvariable=var)
        val.grid(row=row, column=2, sticky="w", padx=4)

    def _build_plot(self):
        self.fig, self.ax = plt.subplots(figsize=(7.2, 4.8), dpi=100)
        self.ax.set_xlabel("Distance (m)")
        self.ax.set_ylabel("Height (m)")
        self.ax.set_title("Trajectory View")
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # keep references to drawable artists & animations
        self.lines = []   # list of Line2D
        self.dots = []    # list of Line2D (markers)
        self.explosions = []  # list of artists
        self.animations = []

    # ----- Utility drawing -----
    def clear_plot(self):
        # Remove artists cleanly
        for ln in self.lines:
            ln.remove()
        for dt in self.dots:
            dt.remove()
        for ex in self.explosions:
            ex.remove()
        self.lines.clear()
        self.dots.clear()
        self.explosions.clear()
        self.animations.clear()

        self.ax.cla()
        self.ax.set_xlabel("Distance (m)")
        self.ax.set_ylabel("Height (m)")
        title = "Bar Chart View" if self.current_mode.get() == "Bar Chart" else "Trajectory View"
        self.ax.set_title(title)
        self.ax.grid(True)
        self.canvas.draw_idle()

    def clear_history(self):
        self.history.clear()
        self.last_race.clear()
        self.stats_text.delete("1.0", "end")
        self.redraw()

    def log_live(self, text):
        self.stats_text.insert("end", text + "\n")
        self.stats_text.see("end")

    # ----- Bar chart rendering -----
    def draw_bar_chart(self, series, title="Results"):
        """
        series: list of dicts with keys {'name','stats','color'}
        Shows two bars per series (max_height, distance).
        """
        self.clear_plot()
        if not series:
            self.ax.text(0.5, 0.5, "No data yet", ha="center", va="center", transform=self.ax.transAxes)
            self.canvas.draw_idle()
            return

        labels = [s['name'] for s in series]
        heights = [s['stats']['max_height'] for s in series]
        dists = [s['stats']['distance'] for s in series]
        xs = range(len(series))

        # bar positions
        width = 0.35
        xs1 = [x - width/2 for x in xs]
        xs2 = [x + width/2 for x in xs]

        # Match colors where possible
        bar1 = self.ax.bar(xs1, heights, width=width, label="Max Height (m)")
        bar2 = self.ax.bar(xs2, dists, width=width, label="Distance (m)")

        self.ax.set_xticks(list(xs))
        self.ax.set_xticklabels(labels, rotation=0)
        self.ax.set_ylabel("Meters")
        self.ax.set_title(title)
        self.ax.legend()
        self.canvas.draw_idle()

    # ----- Redraw depending on mode -----
    def redraw(self):
        mode = self.current_mode.get()
        if mode == "Bar Chart":
            # Prefer last race for comparison; fallback to history last few runs
            data = self.last_race if self.last_race else self.history[-4:]
            self.draw_bar_chart(data, title="Comparison")
        else:
            # Just replot ghosts (history) without animation
            self.clear_plot()
            for item in self.history[-5:]:  # limit clutter
                xs = [p[0] for p in item['positions']]
                ys = [p[1] for p in item['positions']]
                (line,) = self.ax.plot(xs, ys, color=item['color'], alpha=0.35, linewidth=2, label=item['name'])
                self.lines.append(line)
            if self.lines:
                self.ax.legend(loc="upper right")
            self.canvas.draw_idle()

    # ----- Simulation buttons -----
    def launch_one(self):
        if self.current_mode.get() == "Bar Chart":
            self.current_mode.set("Trajectory")

        self.clear_plot()
        fin = self.fin_var.get()
        cd = FINS[fin]
        mass = float(self.mass_var.get())
        angle = float(self.angle_var.get())
        speed = float(self.speed_var.get())
        wind_on = bool(self.wind_var.get())

        # Simulate once, store results
        dt = 0.02
        positions, stats = simulate_flight(mass, cd, angle, speed, dt=dt, wind_on=wind_on)
        color = FIN_COLORS.get(fin, "black")
        name = f"{fin} ({angle:.0f}°, {speed:.0f}m/s)"
        self.history.append(dict(name=name, positions=positions, stats=stats, color=color))

        if not positions:
            messagebox.showwarning("No Flight", "The rocket did not gain altitude. Try higher speed/angle.")
            return

        xs = [p[0] for p in positions]
        ys = [p[1] for p in positions]

        (line,) = self.ax.plot([], [], '-', linewidth=2.5, color=color, label=name)
        (dot,) = self.ax.plot([], [], 'o', color=color, markersize=6)
        self.lines.append(line)
        self.dots.append(dot)
        self.ax.legend(loc="upper right")

        # Live stats
        self.stats_text.delete("1.0", "end")
        self.log_live(f"Single Launch → {name}")
        self.log_live(f"Mass={mass:.2f} kg | Cd≈{cd:.2f} | Wind={'ON' if wind_on else 'OFF'}")

        # Explosion marker (draw when it lands)
        def draw_explosion():
            ex = self.ax.scatter([xs[-1]], [ys[-1] if ys[-1] > 0 else 0], s=80, marker='x', c='red')
            self.explosions.append(ex)

        def init():
            line.set_data([], [])
            dot.set_data([], [])
            return line, dot

        def update(frame):
            line.set_data(xs[:frame+1], ys[:frame+1])
            dot.set_data([xs[frame]], [ys[frame]])
            # live update
            t_now = frame * dt
            self.stats_text.delete("1.0", "end")
            self.log_live(f"t = {t_now:.2f} s | x = {xs[frame]:.1f} m | y = {ys[frame]:.1f} m")

            # axes autoscale a bit
            self.ax.set_xlim(0, max(xs) * 1.1)
            self.ax.set_ylim(0, max(ys) * 1.2)

            return line, dot

        ani = animation.FuncAnimation(
            self.fig, update, frames=len(xs),
            init_func=init, blit=True, interval=max(10, int(dt*1000)), repeat=False
        )
        self.animations.append(ani)
        self.canvas.draw_idle()

        # end summary when animation finishes
        def on_finish():
            draw_explosion()
            self.canvas.draw_idle()
            messagebox.showinfo(
                "Flight Summary",
                f"{name}\n"
                f"Max Height: {stats['max_height']:.1f} m\n"
                f"Distance:   {stats['distance']:.1f} m\n"
                f"Time:       {stats['time_of_flight']:.2f} s"
            )
        # schedule a finish callback after animation time
        total_ms = max(1, int(len(xs) * dt * 1000))
        self.root.after(total_ms + 50, on_finish)

    def race_mode(self):
        if self.current_mode.get() == "Bar Chart":
            self.current_mode.set("Trajectory")
        self.clear_plot()
        self.stats_text.delete("1.0", "end")

        mass = float(self.mass_var.get())
        angle = float(self.angle_var.get())
        speed = float(self.speed_var.get())
        wind_on = bool(self.wind_var.get())
        dt = 0.02

        # Sim all fins
        race_results = []
        max_len = 0
        trajectories = {}
        for fin, cd in FINS.items():
            pos, st = simulate_flight(mass, cd, angle, speed, dt=dt, wind_on=wind_on)
            color = FIN_COLORS.get(fin, "black")
            name = f"{fin}"
            race_results.append(dict(name=name, positions=pos, stats=st, color=color))
            trajectories[fin] = pos
            if len(pos) > max_len:
                max_len = len(pos)

        # Save race for bar chart view and replay
        self.last_race = race_results
        self.history.extend(race_results)

        # Pre-create artists
        lines = {}
        dots = {}
        for item in race_results:
            xs = [p[0] for p in item['positions']]
            ys = [p[1] for p in item['positions']]
            (ln,) = self.ax.plot([], [], '-', linewidth=2.2, color=item['color'], label=item['name'])
            (dtm,) = self.ax.plot([], [], 'o', color=item['color'], markersize=5)
            lines[item['name']] = ln
            dots[item['name']] = dtm
            self.lines.append(ln)
            self.dots.append(dtm)

        self.ax.legend(loc="upper right")

        # victory banner (added at end)
        banner_text = [None]  # use list for mutable closed-over var

        def init():
            for name in lines:
                lines[name].set_data([], [])
                dots[name].set_data([], [])
            return list(lines.values()) + list(dots.values())

        def update(frame):
            # update each rocket to 'frame'
            for item in race_results:
                name = item['name']
                pos = item['positions']
                if frame < len(pos):
                    xs = [p[0] for p in pos[:frame+1]]
                    ys = [p[1] for p in pos[:frame+1]]
                    lines[name].set_data(xs, ys)
                    dots[name].set_data([xs[-1]], [ys[-1]])

            # Bounds
            all_x = [p[0] for it in race_results for p in it['positions'][:frame+1]] or [1]
            all_y = [p[1] for it in race_results for p in it['positions'][:frame+1]] or [1]
            self.ax.set_xlim(0, max(all_x) * 1.1)
            self.ax.set_ylim(0, max(all_y) * 1.2)

            # live text: show best current x
            lead = max(race_results, key=lambda it: (it['positions'][frame][0] if frame < len(it['positions']) else -1))
            lead_x = lead['positions'][frame][0] if frame < len(lead['positions']) else 0.0
            self.stats_text.delete("1.0", "end")
            self.log_live(f"Race: angle {angle:.0f}°, speed {speed:.0f} m/s | Wind={'ON' if wind_on else 'OFF'}")
            self.log_live(f"Leader: {lead['name']} at x = {lead_x:.1f} m")

            return list(lines.values()) + list(dots.values())

        ani = animation.FuncAnimation(
            self.fig, update, frames=max_len,
            init_func=init, blit=True, interval=max(12, int(dt*1000)), repeat=False
        )
        self.animations.append(ani)
        self.canvas.draw_idle()

        # After finish: explosions + banner
        total_ms = max(1, int(max_len * dt * 1000))

        def on_finish():
            # explosions at landing
            for it in race_results:
                if it['positions']:
                    x_last, y_last = it['positions'][-1]
                    ex = self.ax.scatter([x_last], [max(0, y_last)], s=70, marker='x', c='red')
                    self.explosions.append(ex)

            # winner by farthest distance
            winner = max(race_results, key=lambda it: it['stats']['distance'])
            if banner_text[0] is None:
                banner_text[0] = self.ax.text(
                    0.5, 0.95, f"{winner['name'].title()} Wins 🏆",
                    transform=self.ax.transAxes, ha="center", va="top",
                    fontsize=14, weight="bold", color="black",
                    bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.8)
                )
            self.canvas.draw_idle()

            # End summary popup
            summary = "\n".join(
                f"{it['name'].title()}: H={it['stats']['max_height']:.1f}m, D={it['stats']['distance']:.1f}m, T={it['stats']['time_of_flight']:.2f}s"
                for it in race_results
            )
            messagebox.showinfo("Race Summary", summary)

        self.root.after(total_ms + 50, on_finish)

    def replay_ghosts(self):
        # Show last few trajectories as faint ghosts
        self.current_mode.set("Trajectory")
        self.clear_plot()
        ghosts = self.history[-6:]
        if not ghosts:
            messagebox.showinfo("Replay", "No previous launches to replay yet.")
            return
        for item in ghosts:
            xs = [p[0] for p in item['positions']]
            ys = [p[1] for p in item['positions']]
            (ln,) = self.ax.plot(xs, ys, color=item['color'], alpha=0.35, linewidth=2.2, label=item['name'])
            self.lines.append(ln)
        self.ax.legend(loc="upper right")
        self.canvas.draw_idle()


def main():
    root = tk.Tk()
    # Use a nicer default theme where available
    try:
        root.call("source", "sun-valley.tcl")  # optional if you have a ttk theme file
        root.call("set_theme", "light")
    except Exception:
        pass
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
