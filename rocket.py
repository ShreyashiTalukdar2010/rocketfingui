import math

class Rocket:
    def __init__(self, mass, fin_shape, drag_coeff, lift_coeff):
        self.mass = mass
        self.fin_shape = fin_shape
        self.drag_coeff = drag_coeff
        self.lift_coeff = lift_coeff
        self.x = 0
        self.y = 0
        self.vx = 20  # initial horizontal speed
        self.vy = 40  # initial vertical speed
        self.g = 9.81
        self.angle = math.radians(90)  # straight up

    def update(self, dt):
        # Update forces
        speed = math.hypot(self.vx, self.vy)
        if speed == 0:
            return  # rocket’s not moving yet

        # Direction unit vector
        direction_x = self.vx / speed
        direction_y = self.vy / speed

        # Drag force
        drag_force = 0.5 * self.drag_coeff * speed**2
        fx_drag = -drag_force * direction_x
        fy_drag = -drag_force * direction_y

        # Gravity force
        fy_gravity = -self.mass * self.g

        #LIFT force 
        lift_force = 0.5 * self.lift_coeff * speed**2
        fx_lift = -lift_force * direction_y
        fy_lift = lift_force * direction_x

        # Total forces
        fx = fx_drag + fx_lift
        fy = fy_drag + fy_gravity + fy_lift

        # Acceleration
        ax = fx / self.mass
        ay = fy / self.mass

        # Velocity update
        self.vx += ax * dt
        self.vy += ay * dt

        # Position update
        self.x += self.vx * dt
        self.y += self.vy * dt


