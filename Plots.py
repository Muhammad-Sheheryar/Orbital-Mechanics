import json
import os                         # <-- added for folder handling
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# 0.  Common settings
# -------------------------------------------------------------------
# Ensure the "graphs" folder exists before saving anything
os.makedirs("graphs", exist_ok=True)

methods = ["euler_ns", "euler_s", "leapfrog", "rk4"]
labels  = ["Explicit Euler", "Semi‑implicit Euler", "Leapfrog", "RK4"]
colors  = ["red", "green", "blue", "orange"]
plt.rcParams.update({'font.size': 10, 'lines.linewidth': 1.5})

# -------------------------------------------------------------------
# 1.  SEPARATE trajectory sub‑plots (2×2 grid)
# -------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 10))
fig.suptitle("Earth's Orbital Trajectory – One Integrator per Panel", fontsize=14)

for ax, method, label, color in zip(axes.flat, methods, labels, colors):
    with open(f"orbit_data_{method}.json") as f:
        data = json.load(f)

    x = [d["earth_x"] for d in data]
    y = [d["earth_y"] for d in data]
    ax.plot(x, y, color=color, alpha=0.8, linewidth=1.2)
    ax.scatter(0, 0, marker='*', s=150, color='yellow',
               edgecolors='black', zorder=10)
    ax.set_title(label, fontsize=11)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("graphs/orbital_trajectories_separate.pdf", bbox_inches="tight")   # PDF
plt.show()

# -------------------------------------------------------------------
# 2.  Energy error  (log scale, all methods together)
# -------------------------------------------------------------------
plt.figure(figsize=(10, 5.5))
for method, label, color in zip(methods, labels, colors):
    with open(f"orbit_data_{method}.json") as f:
        data = json.load(f)
    t = [d["time"] for d in data]
    e_err = [d["energy_error"] for d in data]
    plt.plot(t, e_err, color=color, alpha=0.85, label=label)

plt.xlabel("Time (s)")
plt.ylabel("Relative Energy Error  |ΔE / E₀|")
plt.title("Total Energy Conservation Error")
plt.legend(fontsize=9)
plt.grid(True, alpha=0.3, which='both')
plt.yscale("log")
plt.tight_layout()
plt.savefig("graphs/energy_error_clean.pdf", bbox_inches="tight")             # PDF
plt.show()

# -------------------------------------------------------------------
# 3.  Angular momentum error
# -------------------------------------------------------------------
plt.figure(figsize=(10, 5.5))
for method, label, color in zip(methods, labels, colors):
    with open(f"orbit_data_{method}.json") as f:
        data = json.load(f)
    t = [d["time"] for d in data]
    l_err = [d["angular_momentum_relative_error"] for d in data]
    plt.plot(t, l_err, color=color, alpha=0.85, label=label)

plt.xlabel("Time (s)")
plt.ylabel("Relative Angular Momentum Error  |ΔLz / Lz₀|")
plt.title("Total Angular Momentum Conservation Error")
plt.legend(fontsize=9)
plt.grid(True, alpha=0.3, which='both')
plt.yscale("log")
plt.tight_layout()
plt.savefig("graphs/angular_momentum_error_clean.pdf", bbox_inches="tight")   # PDF
plt.show()

# -------------------------------------------------------------------
# 4.  Orbital distance stability
# -------------------------------------------------------------------
plt.figure(figsize=(10, 5.5))
for method, label, color in zip(methods, labels, colors):
    with open(f"orbit_data_{method}.json") as f:
        data = json.load(f)
    t = [d["time"] for d in data]
    r = [d["distance_from_star"] for d in data]
    plt.plot(t, r, color=color, alpha=0.85, label=label)

plt.xlabel("Time (s)")
plt.ylabel("Sun–Earth distance (m)")
plt.title("Orbital Distance Stability")
plt.legend(fontsize=9)
plt.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig("graphs/distance_stability_clean.pdf", bbox_inches="tight")       # PDF
plt.show()


# -------------------------------------------------------------------
# 5.  PHASE SPACE PORTRAITS (x vs v_x) - INDIVIDUAL SCALING
# -------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 10))
fig.suptitle("Phase Space Portraits: x vs v_x", fontsize=14)

for ax, method, label, color in zip(axes.flat, methods, labels, colors):
    with open(f"orbit_data_{method}.json") as f:
        data = json.load(f)

    x = [d["earth_x"] for d in data]
    vx = [d["earth_vx"] for d in data]

    ax.plot(x, vx, color=color, alpha=0.7, linewidth=0.8)
    ax.set_title(label, fontsize=11)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("v_x (m/s)")
    ax.grid(True, alpha=0.3)
    ax.set_aspect('auto')  # Allow the ellipse to look like an ellipse, not a circle

plt.tight_layout()
plt.savefig("graphs/phase_space_x_vx_individual.pdf", bbox_inches="tight")
plt.show()