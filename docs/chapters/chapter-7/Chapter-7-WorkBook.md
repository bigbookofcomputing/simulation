## ⚙️ Chapter 7: Physics III: Molecular Dynamics (MD) (Workbook)

The goal of this chapter is to introduce the simulation of **real-time motion** by integrating Newton's equations, contrasting this dynamic approach with the equilibrium sampling of Monte Carlo (MC).

| Section | Topic Summary |
| :--- | :--- |
| **7.1** | Chapter Opener: From Sampling to Dynamics |
| **7.2** | The Velocity–Verlet Algorithm |
| **7.3** | Periodic Boundary Conditions and Neighbor Lists |
| **7.4** | Thermostats and Ensembles (NVE, NVT, NPT) |
| **7.5** | Computing Observables (Energy, Pressure, Diffusion, and Correlation Functions) |

***

### 7.1 From Sampling to Dynamics

> **Summary:** Molecular Dynamics (MD) simulates **continuous time evolution** by integrating Newton's equation, $m_i \frac{d^2 \mathbf{r}_i}{dt^2} = \mathbf{F}_i$. In MD, the **force** $\mathbf{F} = -\nabla E$ replaces the Boltzmann **probability** $e^{-\beta E}$ as the fundamental driving mechanism.

#### Section Detail

MD is a "different worldview" from Monte Carlo. While MC explores statistical equilibrium through random jumps, MD simulates deterministic, time-reversible trajectories driven by the force field of the system's potential energy landscape. The core integrator for MD is the **Velocity-Verlet algorithm**, which is **symplectic** (conserves phase space geometry) and ensures long-term energy stability.

#### Quiz Questions

**1. What is the fundamental driving mechanism for particle movement in a Molecular Dynamics simulation?**

* **A.** Random acceptance/rejection rules based on Boltzmann probability.
* **B.** **The force $\mathbf{F}_i$, derived from the gradient of the potential energy $E(\mathbf{r})$**. (**Correct**)
* **C.** The autocorrelation function.
* **D.** The time step $\Delta t$.

**2. The primary reason the Velocity–Verlet algorithm is preferred over the simpler Euler method for long-term MD simulations is that Velocity–Verlet is:**

* **A.** Faster to compute.
* **B.** Easier to code.
* **C.** **Symplectic and time-reversible, leading to excellent long-term energy conservation**. (**Correct**)
* **D.** A third-order accurate integrator.

---

#### Interview-Style Question

**Question:** Explain the core conceptual difference between the output of a standard Metropolis Monte Carlo simulation and a Molecular Dynamics simulation, even when both model the same system (e.g., liquid Argon).

**Answer Strategy:**
* **MC Output:** Provides a set of configurations weighted by $e^{-\beta E}$ (the **equilibrium ensemble**). The output is sufficient for calculating thermodynamic averages (like $\langle E \rangle$ or specific heat), but it has **no concept of time or dynamics**.
* **MD Output:** Provides a **time-dependent trajectory** of positions and velocities $(\mathbf{r}(t), \mathbf{v}(t))$. This output allows calculation of **dynamic and transport properties** (like diffusion coefficients and correlation functions) that are completely inaccessible to MC.

---

***

### 7.2 The Velocity–Verlet Algorithm

> **Summary:** The **Velocity–Verlet algorithm** is a **second-order accurate** integrator that updates position using the current state and velocity using the average of the current and future forces. It ensures long-term stability by conserving the geometry of phase space (symplectic property).

#### Section Detail

Velocity–Verlet discretizes Newton's equations by splitting the calculation into three sequential steps: position update, new force evaluation, and velocity update. The second-order accuracy ($\mathcal{O}(\Delta t^2)$) ensures that errors scale well with the chosen time step $\Delta t$. The choice of $\Delta t$ is critical and must be small enough to resolve the fastest oscillations in the system (e.g., bond stretches).

#### Quiz Questions

**1. The Velocity–Verlet algorithm uses which two physical quantities at the next time step $t+\Delta t$ to update the velocity $\mathbf{v}_i(t+\Delta t)$?**

* **A.** Position $\mathbf{r}(t+\Delta t)$ and potential energy $U(t+\Delta t)$.
* **B.** The total energy $E_{\text{tot}}(t)$ and the kinetic energy $K(t)$.
* **C.** **The current force $\mathbf{F}(t)$ and the new force $\mathbf{F}(t+\Delta t)$**. (**Correct**)
* **D.** The pressure $P(t)$ and the temperature $T(t)$.

**2. A large MD time step $\Delta t$ that fails to resolve the fastest oscillations in the system primarily leads to:**

* **A.** Incorrect ensemble sampling.
* **B.** **Numerical instability and energy drift**. (**Correct**)
* **C.** High autocorrelation times.
* **D.** Inaccurate pressure calculation.

---

#### Interview-Style Question

**Question:** Briefly explain the "Kick-Drift-Kick" conceptual analogy for the Velocity–Verlet algorithm.

**Answer Strategy:** The Velocity–Verlet algorithm can be viewed as splitting the movement into sequential steps:
1.  **Kick (Half-Step Velocity):** The velocity is advanced by a half-step using the initial acceleration (force).
2.  **Drift (Full-Step Position):** The position is advanced by a full step using this half-step velocity (the drift).
3.  **Kick (Final Velocity):** The force is recalculated at the new position, and the velocity is given its final half-step kick using the average acceleration of the two steps.

---

***

### 7.3 Periodic Boundary Conditions and Neighbor Lists

> **Summary:** **Periodic Boundary Conditions (PBCs)** and the **Minimum Image Convention** are used to emulate an infinite system using a small, finite simulation box. To achieve $\mathcal{O}(N)$ scaling for short-range forces, **Neighbor Lists** are employed to avoid the computationally prohibitive $\mathcal{O}(N^2)$ summation of all pairwise interactions.

#### Section Detail

PBCs eliminate unphysical surface effects by making the simulation box topologically equivalent to a torus. The Minimum Image Convention ensures that each particle interacts only with the nearest periodic image of every other particle. Neighbor lists store pairs within a cutoff radius $r_c$ plus a "skin" buffer $\delta$, and are updated only periodically, dramatically reducing computational cost.

#### Quiz Questions

**1. The **Minimum Image Convention** is used in MD with PBCs to ensure that:**

* **A.** All particles remain in the center of the box.
* **B.** **Each particle interacts only with the nearest periodic image of every other particle**. (**Correct**)
* **C.** The potential energy is always zero.
* **D.** The temperature is constant.

**2. The primary reason for using **Neighbor Lists** in a short-range MD simulation is to reduce the computational complexity of the force calculation from $\mathcal{O}(N^2)$ to approximately:**

* **A.** $\mathcal{O}(\log N)$.
* **B.** $\mathcal{O}(N^3)$.
* **C.** **$\mathcal{O}(N)$**. (**Correct**)
* **D.** $\mathcal{O}(\Delta t)$.

---

#### Interview-Style Question

**Question:** If you are simulating a system using a potential energy function that contains a long-range Coulombic ($1/r$) term, would you still use the simple Neighbor List optimization? Why or why not?

**Answer Strategy:** **No**, the simple neighbor list optimization would be ineffective. The $1/r$ Coulombic potential decays too slowly with distance. Since the interaction cannot be cut off at a finite radius $r_c$ without introducing large errors, every particle still needs to interact with every other particle (and all their images). For long-range forces, specialized methods like the **Ewald summation** or **Particle Mesh Ewald (PME)**, which handle the summation over infinite images, must be used instead of simple cutoffs.

---

***

### 7.4 Thermostats and Ensembles (NVE, NVT, NPT)

> **Summary:** MD allows simulation of different thermodynamic **ensembles**. The simplest is the **NVE** (Microcanonical) ensemble, which conserves energy. To control temperature, a **thermostat** is added (creating the **NVT** ensemble), with the **Nosé–Hoover thermostat** being preferred for generating the statistically correct canonical distribution. The **NPT** ensemble adds a **barostat** to control pressure.

#### Section Detail

* **NVE:** Pure Velocity-Verlet, ideal for testing energy conservation.
* **NVT:** Requires a thermostat to adjust velocities and keep the instantaneous temperature $T_{\text{inst}}$ close to the target $T_0$. The Berendsen thermostat is simple (velocity rescaling), while Nosé–Hoover introduces auxiliary dynamics ($\xi$) to accurately sample the canonical ensemble.
* **NPT:** Allows the simulation box volume $V$ to fluctuate to maintain a target pressure $P_0$.

#### Quiz Questions

**1. The goal of the **Nosé–Hoover thermostat** is to ensure the MD simulation correctly samples which thermodynamic ensemble?**

* **A.** The Microcanonical (NVE) ensemble.
* **B.** The Isobaric-Isothermal (NPT) ensemble.
* **C.** **The Canonical (NVT) ensemble**. (**Correct**)
* **D.** The Grand Canonical ($\mu V T$) ensemble.

**2. To simulate a system at constant temperature ($T$) and constant pressure ($P$), which type of simulation must be run?**

* **A.** Monte Carlo simulation.
* **B.** NVE simulation.
* **C.** **NPT simulation (Isothermal–Isobaric) using both a thermostat and a barostat**. (**Correct**)
* **D.** NVT simulation.

---

#### Interview-Style Question

**Question:** Why is the simple Berendsen thermostat often only used for the **equilibration** phase of an MD simulation, and not for the final **production** (measurement) phase?

**Answer Strategy:** The Berendsen thermostat achieves temperature control by non-physically rescaling velocities based on the difference between the instantaneous temperature and the target temperature. While this method is robust for quickly bringing the system to the target $T$ (equilibration), it does **not generate the statistically correct canonical ensemble**. Specifically, it suppresses energy fluctuations, which means observables calculated during the production phase (e.g., specific heat, which depends on energy fluctuations) will be inaccurate. For production, a method like Nosé–Hoover or Langevin is required.

---

***

### 7.5 Computing Observables

> **Summary:** MD extracts physics from trajectories by computing **time-averaged observables**. **Pressure** requires the **Virial Theorem** (mixing kinetic and inter-particle force terms). **Transport properties** are calculated from time correlation functions, such as the **Mean-Squared Displacement (MSD)** or the **Velocity Autocorrelation Function (VACF)**, both of which yield the **Diffusion Coefficient**.

#### Section Detail

The MSD measures the average distance a particle travels from its origin, $\text{MSD}(t) = \langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle$, which is linear in time for diffusive systems: $D = \lim_{t \to \infty} \frac{1}{6t} \text{MSD}(t)$. The VACF, $C_v(t) = \frac{\langle \mathbf{v}(0) \cdot \mathbf{v}(t) \rangle}{\langle \mathbf{v}(0)^2 \rangle}$, captures the system's memory of motion. The total energy $E_{\text{tot}} = K+U$ serves as the primary diagnostic.

#### Quiz Questions

**1. Which theorem is used in MD to calculate the system's pressure $P$, by including a contribution from the interparticle forces $\mathbf{F}_{ij}$?**

* **A.** The Fluctuation-Dissipation Theorem.
* **B.** The Equipartition Theorem.
* **C.** **The Virial Theorem**. (**Correct**)
* **D.** The Intermediate Value Theorem.

**2. Which two time-dependent functions are used to calculate the Diffusion Coefficient ($D$)?**

* **A.** Total Energy $E(t)$ and Pressure $P(t)$.
* **B.** The Pressure-Volume product $PV$ and $N k_B T$.
* **C.** **The Mean-Squared Displacement (MSD) and the Velocity Autocorrelation Function (VACF)**. (**Correct**)
* **D.** The potential energy $U(t)$ and the time step $\Delta t$.

---

#### Interview-Style Question

**Question:** A physicist simulates a liquid and measures the Velocity Autocorrelation Function, $C_v(t)$. They notice that $C_v(t)$ initially drops quickly but then becomes slightly negative before decaying to zero. Explain the physical origin of this negative correlation.

**Answer Strategy:** A negative correlation in $C_v(t)$ means the particle, after a short time $\tau$, is more likely to be moving **in the opposite direction** ($\mathbf{v}(0) \cdot \mathbf{v}(\tau) < 0$).
* This is characteristic of a **liquid** or dense fluid.
* The negative value is caused by **caging effects**: a central particle, initially moving at $\mathbf{v}(0)$, collides with its dense, surrounding shell of neighbors. The particle "bounces" off the surrounding cage, reversing its initial velocity vector, causing the instantaneous velocity to correlate negatively with the initial velocity.

---

***

## 💡 Hands-On Simulation Projects (Chapter Conclusion) 🛠️

These projects are designed to implement the core MD pipeline and compute essential equilibrium and dynamic properties.

### Project 1: Implementing the Velocity–Verlet Integrator (The Engine)

* **Goal:** Implement the core Velocity–Verlet algorithm for a single particle in a 1D quadratic potential $U(r) = \frac{1}{2} k r^2$ (Harmonic Oscillator).
* **Setup:** Use $m=1.0, k=1.0, \Delta t=0.01$. Initial conditions $r_0=1.0, v_0=0.0$.
* **Steps:**
    1.  Define the force function $\mathbf{F}(r) = -k r$.
    2.  Implement the full three-step Velocity–Verlet update loop (position $\to$ new force $\to$ velocity).
    3.  Run for $5000$ steps and record the total energy $E_{\text{tot}} = K + U$ at each step.
* ***Goal***: Plot $E_{\text{tot}}(t)$ and show that it remains constant (with only small numerical fluctuations), confirming the symplectic stability of the integrator.

### Project 2: MD with Periodic Boundaries and Collision

* **Goal:** Extend the simulation to 2D with PBC and non-trivial interactions (conceptual Lennard-Jones).
* **Setup:** Place $N=4$ particles in a square box of side $L=10.0$. Set up initial random positions and zero velocity.
* **Steps:**
    1.  Implement the `minimum_image` function for calculating the shortest distance between two points under PBC.
    2.  Define a conceptual repulsive force: $\mathbf{F}_{ij} \propto 1/r^7$ for $r < 1.0$, and zero otherwise.
    3.  Run the Velocity–Verlet loop, applying the **PBC wrapping** to the positions $\mathbf{r}_i$ after every full step.
* ***Goal***: Demonstrate particle movement and successful wrapping when particles cross the box boundaries.

### Project 3: Computing the Diffusion Coefficient ($D$)

* **Goal:** Calculate the diffusion coefficient by measuring the Mean-Squared Displacement (MSD).
* **Setup:** Simulate a system of particles ($N \gg 1$) at a high temperature (liquid/gas state). Record the positions $\mathbf{r}(t)$ at regular intervals over a long trajectory.
* **Steps:**
    1.  Calculate the MSD over the time trajectory: $\text{MSD}(\tau) = \langle |\mathbf{r}(t+\tau) - \mathbf{r}(t)|^2 \rangle$.
    2.  Plot $\text{MSD}(\tau)$ versus $\tau$.
    3.  Fit the long-time, linear regime of the MSD curve to a straight line: $\text{MSD}(\tau) = 6 D \tau + C$.
    4.  Extract the slope and compute the diffusion coefficient $D$.
* ***Goal***: Confirm the expected linear growth of MSD in a diffusive system and obtain a quantitative transport property.

### Project 4: Implementing the Berendsen Thermostat (NVT)

* **Goal:** Modify the NVE integrator to simulate a canonical (NVT) ensemble by controlling temperature.
* **Setup:** Use the same Harmonic Oscillator (Project 1) but initialize with a higher energy (e.g., $r_0=5.0, v_0=0.0$). Set a target temperature $T_0=1.0$ and relaxation time $\tau_T=1.0$.
* **Steps:**
    1.  At each step, calculate the instantaneous temperature $T_{\text{inst}}$ from the kinetic energy $K$.
    2.  Calculate the Berendsen velocity scaling factor $\lambda$.
    3.  Apply the scaling: $\mathbf{v} \leftarrow \lambda \mathbf{v}$ immediately before the next position update.
* ***Goal***: Plot the instantaneous temperature $T_{\text{inst}}(t)$ and show that it smoothly and quickly relaxes from the high initial temperature to the target temperature $T_0$, demonstrating the successful implementation of the NVT control.


