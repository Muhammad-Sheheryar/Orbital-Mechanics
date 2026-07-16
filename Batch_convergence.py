import os
import math
import numpy as np
import matplotlib.pyplot as plt

# Import your simulation classes – adjust if your main file has a different name.
from simulation import OrbitalSimulation, create_earth_sun_system

# ----------------------------------------------------------------------
# Extract constants DIRECTLY from the simulator to guarantee exact match
# ----------------------------------------------------------------------
# Create a dummy simulation just to get the initial bodies
dummy_sim = OrbitalSimulation("rk4", 1.0)  
sun, earth = create_earth_sun_system()

# Masses
M_sun = sun.mass
M_earth = earth.mass
G = OrbitalSimulation.G

# Initial relative position (Earth relative to Sun)
r0 = math.sqrt((earth.x - sun.x)**2 + (earth.y - sun.y)**2)

# Exact two-body angular frequency (using reduced mass)
omega = np.sqrt(G * (M_sun + M_earth) / r0**3)

# Exact orbital period
T = 2 * np.pi / omega

# Integration time
TARGET_ORBITS = 1.5
t_final = TARGET_ORBITS * T

print(f"Extracted constants:")
print(f"  M_sun   = {M_sun:.4e} kg")
print(f"  M_earth = {M_earth:.4e} kg")
print(f"  r0      = {r0:.4e} m")
print(f"  omega   = {omega:.4e} rad/s")
print(f"  T       = {T:.4e} s  ({T/86400:.2f} days)")
print(f"  t_final = {t_final:.4e} s  ({TARGET_ORBITS} orbits)")

# ----------------------------------------------------------------------
# Analytical reference solution (exact for the two-body circular orbit)
# ----------------------------------------------------------------------
def get_reference_state():
    """
    Returns the exact relative position (x_rel, y_rel) of the Earth
    with respect to the Sun at t = t_final.
    
    """
    x_rel_true = r0 * np.cos(omega * t_final)
    y_rel_true = r0 * np.sin(omega * t_final)
    return x_rel_true, y_rel_true

# ----------------------------------------------------------------------
# Helper: run simulation and return final state of BOTH bodies
# ----------------------------------------------------------------------
def run_simulation_orbits(method, dt_frac):
    dt = dt_frac * T
    sim = OrbitalSimulation(method, dt)
    sun, earth = create_earth_sun_system()
    
    # ---- PATCH: Use the exact two-body circular velocity ----
    # The relative distance
    r0 = math.sqrt((earth.x - sun.x)**2 + (earth.y - sun.y)**2)
    # Exact circular velocity for the relative orbit
    v_circ = math.sqrt(G * (M_sun + M_earth) / r0)
    # Set Earth's velocity to be purely tangential (along +y, since Earth starts at +x)
    earth.vx = 0.0
    earth.vy = v_circ
    # ---- End of patch ----
    
    sim.add_body(sun)
    sim.add_body(earth)
    sim.compute_accelerations()
    
    steps = int(t_final / dt)
    for _ in range(steps):
        sim.step()
    
    return sun.x, sun.y, earth.x, earth.y

# ----------------------------------------------------------------------
# Convergence analysis
# ----------------------------------------------------------------------
def convergence_analysis():
    # Get the analytical reference relative position
    x_ref, y_ref = get_reference_state()
    print(f"\nAnalytical reference at t = {TARGET_ORBITS} T:")
    print(f"  x_rel = {x_ref:.3e} m, y_rel = {y_ref:.3e} m")

    # Use the SAME timestep list for all integrators
    dt_list = [0.0005, 0.001, 0.002, 0.005, 0.01]

    methods = {
        "Explicit Euler": "euler_ns",
        "Semi-implicit Euler": "euler_s",
        "Leapfrog": "leapfrog",
        "RK4": "rk4",
    }

    results = {}   # method_name -> list of (dt_frac, relative_error)

    for name, tag in methods.items():
        errors = []
        print(f"\n--- {name} ---")
        
        for dt_frac in dt_list:
            # Run simulation and get both bodies
            sun_x, sun_y, earth_x, earth_y = run_simulation_orbits(tag, dt_frac)
            
            # Relative position from the simulation
            rel_x = earth_x - sun_x
            rel_y = earth_y - sun_y
            
            # Relative error (dimensionless): |r_num - r_true| / r0
            err = math.sqrt((rel_x - x_ref)**2 + (rel_y - y_ref)**2) / r0
            errors.append(err)
            
            print(f"  dt = {dt_frac:.4f} T  ->  rel. error = {err:.3e}")

        results[name] = list(zip(dt_list, errors))

    # ------------------------------------------------------------------
    # Plot: log-log with reference slopes
    # ------------------------------------------------------------------
    plt.figure(figsize=(9, 7))
    
    # Plot the actual data
    for name, data in results.items():
        dt_vals = np.array([d[0] for d in data])
        err_vals = np.array([d[1] for d in data])
        
        log_dt = np.log10(dt_vals)
        log_err = np.log10(err_vals)
        
        # Empirical slope (linear fit in log-log space)
        slope, intercept = np.polyfit(log_dt, log_err, 1)
        
        plt.loglog(dt_vals, err_vals, 'o-', linewidth=2, markersize=8,
                   label=f"{name} (slope = {slope:.2f})")
    
    # Add reference slope lines (1st, 2nd, 4th order)
    # Anchor the reference lines to the Leapfrog data point at dt=0.001
    # (so they overlap nicely with the data)
    dt_ref = np.array([0.0005, 0.01])
    
    # Find the Leapfrog error at dt=0.001 to anchor the 2nd order line
    leapfrog_data = dict(results["Leapfrog"])
    err_at_001 = leapfrog_data[0.001]
    
    # 1st order: slope 1
    err_1st = err_at_001 * (dt_ref / 0.001)**1.0 * 0.5  # scaled to sit near Euler
    # 2nd order: slope 2
    err_2nd = err_at_001 * (dt_ref / 0.001)**2.0
    # 4th order: slope 4, anchored to RK4 at dt=0.01 (since RK4 hits floor at small dt)
    rk4_data = dict(results["RK4"])
    err_at_001_rk4 = rk4_data[0.001]
    err_4th = err_at_001_rk4 * (dt_ref / 0.001)**4.0
    
    plt.loglog(dt_ref, err_1st, '--', color='gray', alpha=0.7, linewidth=1.5, label='1st order (ref)')
    plt.loglog(dt_ref, err_2nd, '--', color='darkgray', alpha=0.7, linewidth=1.5, label='2nd order (ref)')
    plt.loglog(dt_ref, err_4th, '--', color='lightgray', alpha=0.7, linewidth=1.5, label='4th order (ref)')

    plt.xlabel("Timestep fraction $\\Delta t / T$", fontsize=13)
    plt.ylabel(f"Relative position error after {TARGET_ORBITS} orbits", fontsize=13)
    plt.title(f"Convergence test: Earth–Sun system (analytical reference)", fontsize=14)
    plt.grid(True, which="both", linestyle="--", alpha=0.4)
    plt.legend(loc="lower right", fontsize=10)
    plt.tight_layout()

    # ------------------------------------------------------------------
    # Save the figure
    # ------------------------------------------------------------------
    os.makedirs("graphs", exist_ok=True)
    pdf_filename = f"graphs/convergence_plot_{TARGET_ORBITS:.1f}_orbits_analytical.pdf"
    plt.savefig(pdf_filename, bbox_inches="tight")
    print(f"\nPlot saved to {pdf_filename}")
    plt.show()

    # ------------------------------------------------------------------
    # Print summary
    # ------------------------------------------------------------------
    print("\n=== Empirical convergence orders ===")
    print(f"{'Method':<22} {'Slope':<8} {'Expected'}")
    print("-" * 45)
    for name, data in results.items():
        dt_vals = np.array([d[0] for d in data])
        err_vals = np.array([d[1] for d in data])
        slope, _ = np.polyfit(np.log10(dt_vals), np.log10(err_vals), 1)
        
        if name == "Explicit Euler" or name == "Semi-implicit Euler":
            expected = 1
        elif name == "Leapfrog":
            expected = 2
        else:  # RK4
            expected = 4
        
        print(f"{name:<22} {slope:.3f}      {expected}")

# ----------------------------------------------------------------------
if __name__ == "__main__":
    convergence_analysis()