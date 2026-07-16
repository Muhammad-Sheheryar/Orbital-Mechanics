import math
import json
import pygame
import sys

# ==================
# 1. PHYSICS ENGINE 
# ==================

class CelestialBody:
    def __init__(self, x, y, vx, vy, mass, radius, color):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.ax = 0.0
        self.ay = 0.0
        self.mass = mass
        self.radius = radius
        self.color = color
        # Storing history strictly for phase space and orbital plots
        self.history = {"x": [], "y": [], "vx": [], "vy": []}

    def record_state(self):
        self.history["x"].append(self.x)
        self.history["y"].append(self.y)
        self.history["vx"].append(self.vx)
        self.history["vy"].append(self.vy)

class OrbitalSimulation:
    G = 6.67430e-11

    def __init__(self, method, dt):
        self.method = method.lower()
        self.dt = dt
        self.time = 0.0
        self.bodies = []
        self.data_log = []

    def add_body(self, body):
        self.bodies.append(body)

    def compute_accelerations(self):
        for b in self.bodies:
            b.ax, b.ay = 0.0, 0.0

        for i in range(len(self.bodies)):
            for j in range(i + 1, len(self.bodies)):
                b1, b2 = self.bodies[i], self.bodies[j]
                dx, dy = b2.x - b1.x, b2.y - b1.y
                dist_sq = dx**2 + dy**2
                if dist_sq == 0: continue
                
                dist = math.sqrt(dist_sq)
                force = self.G * b1.mass * b2.mass / dist_sq
                
                fx = force * (dx / dist)
                fy = force * (dy / dist)

                b1.ax += fx / b1.mass
                b1.ay += fy / b1.mass
                b2.ax -= fx / b2.mass
                b2.ay -= fy / b2.mass

    def step(self):
        if self.method == "euler_ns":
            self._step_euler_NONSYMPLECTIC()
        elif self.method == "euler_s":
            self._step_euler_SYMPLECTIC()
        elif self.method == "leapfrog":
            self._step_leapfrog()
        elif self.method == "rk4":
            self._step_rk4()
        
        self.time += self.dt
        for b in self.bodies:
            b.record_state()

    def _step_euler_NONSYMPLECTIC(self):
        self.compute_accelerations()
        for b in self.bodies:
            b.x += b.vx * self.dt
            b.y += b.vy * self.dt
            b.vx += b.ax * self.dt
            b.vy += b.ay * self.dt


    def _step_euler_SYMPLECTIC(self):
        self.compute_accelerations()
        for b in self.bodies:
            b.vx += b.ax * self.dt
            b.vy += b.ay * self.dt
            b.x += b.vx * self.dt
            b.y += b.vy * self.dt



    def _step_leapfrog(self):
        # Kick 1
        for b in self.bodies:
            b.vx += 0.5 * b.ax * self.dt
            b.vy += 0.5 * b.ay * self.dt
        # Drift
        for b in self.bodies:
            b.x += b.vx * self.dt
            b.y += b.vy * self.dt
        # Update Forces
        self.compute_accelerations()
        # Kick 2
        for b in self.bodies:
            b.vx += 0.5 * b.ax * self.dt
            b.vy += 0.5 * b.ay * self.dt

    def _step_rk4(self):
        """
        Synchronous Runge–Kutta 4th order integration for the whole system.
        All bodies are advanced together through the four trial steps.
        """

        # --- Save the original state of every body ---
        for b in self.bodies:
            b._orig_x, b._orig_y = b.x, b.y
            b._orig_vx, b._orig_vy = b.vx, b.vy

        # --- Stage k1 (slope at the beginning) ---
        self.compute_accelerations()          # uses original positions
        for b in self.bodies:
            b.k1vx = b.ax * self.dt
            b.k1vy = b.ay * self.dt
            b.k1x  = b.vx * self.dt
            b.k1y  = b.vy * self.dt

        # --- Stage k2 (midpoint, using k1) ---
        for b in self.bodies:
            b.x  = b._orig_x  + 0.5 * b.k1x
            b.y  = b._orig_y  + 0.5 * b.k1y
            b.vx = b._orig_vx + 0.5 * b.k1vx
            b.vy = b._orig_vy + 0.5 * b.k1vy
        self.compute_accelerations()          # forces at the common midpoint
        for b in self.bodies:
            b.k2vx = b.ax * self.dt
            b.k2vy = b.ay * self.dt
            b.k2x  = (b._orig_vx + 0.5 * b.k1vx) * self.dt
            b.k2y  = (b._orig_vy + 0.5 * b.k1vy) * self.dt

        # --- Stage k3 (second midpoint, using k2) ---
        for b in self.bodies:
            b.x  = b._orig_x  + 0.5 * b.k2x
            b.y  = b._orig_y  + 0.5 * b.k2y
            b.vx = b._orig_vx + 0.5 * b.k2vx
            b.vy = b._orig_vy + 0.5 * b.k2vy
        self.compute_accelerations()
        for b in self.bodies:
            b.k3vx = b.ax * self.dt
            b.k3vy = b.ay * self.dt
            b.k3x  = (b._orig_vx + 0.5 * b.k2vx) * self.dt
            b.k3y  = (b._orig_vy + 0.5 * b.k2vy) * self.dt

        # --- Stage k4 (endpoint, using k3) ---
        for b in self.bodies:
            b.x  = b._orig_x  + b.k3x
            b.y  = b._orig_y  + b.k3y
            b.vx = b._orig_vx + b.k3vx
            b.vy = b._orig_vy + b.k3vy
        self.compute_accelerations()
        for b in self.bodies:
            b.k4vx = b.ax * self.dt
            b.k4vy = b.ay * self.dt
            b.k4x  = (b._orig_vx + b.k3vx) * self.dt
            b.k4y  = (b._orig_vy + b.k3vy) * self.dt

        # --- Final weighted average (synchronous update) ---
        for b in self.bodies:
            b.x  = b._orig_x  + (b.k1x  + 2*b.k2x  + 2*b.k3x  + b.k4x)  / 6.0
            b.y  = b._orig_y  + (b.k1y  + 2*b.k2y  + 2*b.k3y  + b.k4y)  / 6.0
            b.vx = b._orig_vx + (b.k1vx + 2*b.k2vx + 2*b.k3vx + b.k4vx) / 6.0
            b.vy = b._orig_vy + (b.k1vy + 2*b.k2vy + 2*b.k3vy + b.k4vy) / 6.0

    def calculate_system_energy(self):
            kinetic, potential = 0.0, 0.0
            for i, b in enumerate(self.bodies):
                kinetic += 0.5 * b.mass * (b.vx**2 + b.vy**2)
                for j in range(i + 1, len(self.bodies)):
                    other = self.bodies[j]
                    dist = math.sqrt((other.x - b.x)**2 + (other.y - b.y)**2)
                    potential -= self.G * b.mass * other.mass / dist
            return kinetic + potential
    
    def get_total_energy(self):
        """Returns kinetic + potential energy of the system."""
        return self.calculate_system_energy()

# ==========================================
# 2. INITIALIZATION & DATA LOGGING
# ==========================================

def create_earth_sun_system():
    M_sun, M_earth = 1.98847e30, 5.97237e24
    d = 1.4959826e11
    v_earth = math.sqrt(OrbitalSimulation.G * M_sun / d)

    sun = CelestialBody(0, 0, 0, 0, M_sun, 30, (255, 255, 0))
    earth = CelestialBody(d, 0, 0, v_earth, M_earth, 10, (0, 0, 255))
    return sun, earth
def run_simulation(method, dt_fraction, target_orbits=None,track_power_law=False, total_time=None):
    """Run the math without Pygame. Perfect for batch data generation."""
    period = 365 * 86400
    dt = dt_fraction * period
    
    sim = OrbitalSimulation(method, dt)
    sun, earth = create_earth_sun_system()
    sim.add_body(sun)
    sim.add_body(earth)
    
    sim.compute_accelerations() # Initial forces
    E_initial = sim.calculate_system_energy()



    steps = int((target_orbits * period) / dt)
    L_initial = None   # Will store the first total angular momentum

    for step in range(steps):
        sim.step()
        E_current = sim.calculate_system_energy()
        
        # Calculate total angular momentum
        L_total = (sun.mass * (sun.x * sun.vy - sun.y * sun.vx) +
                earth.mass * (earth.x * earth.vy - earth.y * earth.vx))
        
        # Store initial reference values on first step
        if step == 0:
            L_initial = L_total
            E_initial = E_current   # Also store energy initial here (better to do outside loop)
        
        # Current angular momentum for logging
        L_current = L_total
        
        sim.data_log.append({
            "time": sim.time,
            "energy_error": abs((E_current - E_initial) / E_initial) if E_initial != 0 else 0,
            "angular_momentum": L_current,
            "angular_momentum_relative_error": abs((L_current - L_initial) / L_initial) if L_initial != 0 else 0,
            "distance_from_star": math.sqrt(earth.x**2 + earth.y**2),
            "earth_x": earth.x,
            "earth_y": earth.y,
            "earth_vx": earth.vx,
            "earth_vy": earth.vy
        })

    filename = f"orbit_data_{method}.json"
    with open(filename, "w") as f:
            json.dump(sim.data_log, f, indent=4)
            
    print(f"Success! Saved {len(sim.data_log)} telemetry points to {filename}")

            
     
if __name__ == "__main__":
    run_simulation("euler_ns", 0.01, 100)