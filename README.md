# Orbital-Mechanics: Comparative N-Body Dynamics and Symplectic Geometry

This Python-based orbital simulator models a two-body Newtonian gravitational system to analyze how different numerical integration schemes conserve physical and geometric properties over long timescales. The repository compares four integrators: Explicit Euler, Semi-implicit Euler, Leapfrog (Velocity Verlet), and Runge-Kutta 4 (RK4).

While traditional computational physics prioritizes a solver's formal order of accuracy, this project proves that **symplecticity**—the preservation of phase-space volume—is the critical factor for long-term orbital stability. By tracking an Earth-Sun system over 100 orbits, the simulation reveals that:

* **Non-symplectic methods** (Explicit Euler and RK4) suffer from secular energy drift and angular momentum growth, which eventually degrades the orbital geometry.
* **Symplectic methods** (Semi-implicit Euler and Leapfrog) bound energy errors and prevent secular drift. Leapfrog conserves angular momentum to machine precision ($\sim 10^{-15}$) with minimal computational overhead[cite: 1].

This repository serves as a clean, modular educational toolkit for physicists and developers looking to see Hamiltonian mechanics in action[cite: 1].

---

## 🌌 Data and Visualizations

Since the visualization outputs are generated as high-fidelity vector PDFs, click the links below to view the interactive plots directly[cite: 1]:

* 📈 **[Orbital Trajectories Over 100 Orbits](graphs/orbital_trajectories_separate.pdf):** Compares the dramatic visual consequence of energy inflation in Explicit Euler versus the stable, closed structures of Leapfrog and RK4[cite: 1].
* 🌀 **[Phase-Space Portraits ($x$ vs. $v_x$)](graphs/phase_space_x_vx.pdf):** Direct proof of Hamiltonian structural preservation. Note the closed, invariant curves traced by symplectic methods compared to the outward spiral of Explicit Euler[cite: 1].
* 📊 **[Energy Conservation Analysis](graphs/energy_error_clean.pdf):** Showcases the bounded, oscillating error of symplectic integrators versus the secular, monotonic drift of non-symplectic integrators[cite: 1].
* 📐 **[Angular Momentum Error Plots](graphs/angular_momentum_error_clean.pdf):** Highlights Leapfrog tracking angular momentum conservation down to machine precision limits[cite: 1].
* 📉 **[Integrator Convergence Test](graphs/convergence_plot_1_5_orbits.pdf):** Verifies the empirical scaling order ($\mathcal{O}(\Delta t^n)$) for each solver implementation[cite: 1].
* 📏 **[Distance Stability Metrics](graphs/distance_stability_clean.pdf):** Tracks long-term boundary changes in the semi-major axis over extended timelines[cite: 1].

---

## The Physics & Mathematical Foundations

### 1. Equations of Motion
The simulator models $N$ mutually interacting bodies governed by Newton's Law of Universal Gravitation[cite: 1]:

$$\frac{d^2 \mathbf{r}_i}{dt^2} = \sum_{j \neq i} \frac{G m_j}{\vert{}\mathbf{r}_{ij}\vert{}^3} \mathbf{r}_{ij}$$

Where $\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i$ represents the relative position vector between body $i$ and body $j$[cite: 1]. To solve this numerically, we decouple this second-order system into a set of coupled first-order differential equations[cite: 1]:

$$\frac{d\mathbf{r}_i}{dt} = \mathbf{v}_i, \quad \frac{d\mathbf{v}_i}{dt} = \mathbf{a}_i(\mathbf{r})$$

### 2. Conservation Laws and Geometry
For an exact, isolated two-body system, the total mechanical energy $E$ and total angular momentum vector $\mathbf{L}$ are conserved invariants[cite: 1]:

$$E = \sum_{i} \frac{1}{2} m_i \mathbf{v}_i^2 - \sum_{i < j} \frac{G m_i m_j}{\vert{}\mathbf{r}_{ij}\vert{}}$$

$$\mathbf{L} = \sum_{i} m_i (\mathbf{r}_i \times \mathbf{v}_i)$$

In the $xy$-plane, we track the scalar $z$-component $L_z$[cite: 1]. The conservation of energy stabilizes the semi-major axis $a = -GMm / (2E)$[cite: 1], while the conservation of angular momentum locks in the eccentricity $e$[cite: 1] and satisfies Kepler's second law (constant areal velocity $\frac{dA}{dt} = \frac{L}{2m}$)[cite: 1].

---

## 💻 Implemented Integrators

The mathematical implementation of each numerical engine is detailed below:

### Explicit Euler (1st-Order, Non-Symplectic)
Updates position first, then velocity using the values from the *previous* time step[cite: 1]:

$$\mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{v}_n \Delta t$$

$$\mathbf{v}_{n+1} = \mathbf{v}_n + \mathbf{a}_n \Delta t$$

Because it fails to preserve the symplectic 2-form, it continuously pumps artificial energy into the system, causing the orbit to expand exponentially[cite: 1].

### Semi-Implicit Euler (1st-Order, Symplectic)
Crucially updates the velocity first, then utilizes this *advanced* velocity to update the position[cite: 1]:

$$\mathbf{v}_{n+1} = \mathbf{v}_n + \mathbf{a}_n \Delta t$$

$$\mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{v}_{n+1} \Delta t$$

This simple swap makes the method symplectic[cite: 1]. It preserves phase-space volume, forcing the energy error to remain bounded (oscillating around a "shadow Hamiltonian") rather than drifting secularly[cite: 1].

### Leapfrog / Velocity Verlet (2nd-Order, Symplectic)
A staggered half-step velocity kick, full-step position drift, and second half-step velocity kick[cite: 1]:

$$\mathbf{v}_{n+1/2} = \mathbf{v}_n + \frac{\Delta t}{2} \mathbf{a}_n$$

$$\mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{v}_{n+1/2} \Delta t$$

$$\mathbf{v}_{n+1} = \mathbf{v}_{n+1/2} + \frac{\Delta t}{2} \mathbf{a}_{n+1}$$

Leapfrog is highly computationally efficient (only $1$ force evaluation per step)[cite: 1] and conserves angular momentum down to double-precision machine limits ($\sim 10^{-15}$)[cite: 1].

### Runge-Kutta 4 (4th-Order, Non-Symplectic)
Calculates a weighted average of four trial slopes across the interval $\Delta t$[cite: 1]:

$$\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{\Delta t}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)$$

Where $\mathbf{y} = [\mathbf{r}, \mathbf{v}]^T$ is the full system state vector[cite: 1]. While locally incredibly accurate, its non-symplectic nature leads to linear accumulation of energy error over long integrations[cite: 1].

---

## 📊 Performance and Scaling Analysis

| Integrator | Order of Accuracy | Force Evaluations / Step | Energy Error ($\Delta E / E_0$) | Angular Momentum Error ($\Delta L_z / L_0$) | Long-Term Geometric Fidelity |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Explicit Euler** | $\mathcal{O}(\Delta t)$ | $1$ | Unbounded (Secular Growth)[cite: 1] | $\sim 5 \times 10^{-5}$[cite: 1] | Completely Fails (Spirals Outward)[cite: 1] |
| **Semi-Implicit Euler** | $\mathcal{O}(\Delta t)$ | $1$ | Bounded ($\sim 10^{-3}$)[cite: 1] | $\sim 10^{-14}$[cite: 1] | Bounded with Visible Precession[cite: 1] |
| **Leapfrog** | $\mathcal{O}(\Delta t^2)$ | $1$ | Bounded ($\sim 10^{-5}$)[cite: 1] | $\sim 10^{-15}$ (Machine Limit)[cite: 1] | **Excellent** (No visual precession)[cite: 1] |
| **Runge-Kutta 4** | $\mathcal{O}(\Delta t^4)$ | $4$ | Secular Drift ($< 10^{-7}$)[cite: 1] | $\sim 2 \times 10^{-10}$[cite: 1] | Degrades over long timescales ($10^4+$ orbits)[cite: 1] |

---

## 🛠️ Execution & Pipeline Flow

The codebase separates telemetry simulation from vector visualization to maximize structural efficiency[cite: 1]:

### 1. Compute State Telemetry
To run the primary gravitational engine, execute `simulation.py`[cite: 1]. This calculates the dynamics over 100 orbits and saves raw telemetry state data into individual `.json` matrix arrays within the `data/` folder[cite: 1]:
```bash
python simulation.py
```
## Render Geometric & Conservation Plots
Once the telemetry JSON files are generated, run `plot.py` to parse the telemetry records and render the orbital trajectories, phase-space contours, distance metrics, energy errors, and angular momentum arrays into vector PDFs inside the graphs/ directory[cite: 1]:
```bash
python plot.py
```

## Analyze Power-Law Convergence

To analyze systematic step-size scaling limits, run `Batch_convergence.py`[cite: 1]. This script automatically maps a comprehensive multi-timestep error sweep across all integrated routines, computes empirical error regression coefficients, and self-generates the final log-log power-law convergence graph (graphs/convergence_plot_1_5_orbits.pdf) upon execution[cite: 1]:

```bash
python Batch_convergence.py
``` 
