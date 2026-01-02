##  Chapter 2: Physics I – The Ising Model (Workbook)

The goal of this chapter is to apply the Markov Chain Monte Carlo (MCMC) engine to the **2D Ising Model**, demonstrating how simple local rules give rise to complex emergent phenomena, like the phase transition.

| Section | Topic Summary |
| :--- | :--- |
| **2.1** | Chapter Opener: Emergence and the Grand Simulation |
| **2.2** | The Hamiltonian and the Local Rule |
| **2.3** | The Computational Framework: Periodic Boundary Conditions |
| **2.4** | Implementing the Metropolis Algorithm on the Lattice |
| **2.5** | Analysis I: Macroscopic Observables |
| **2.6** | Analysis II: Equilibration and Autocorrelation |
| **2.7** | Core Application: Locating the Phase Transition |



### 2.1 Emergence and the Grand Simulation

> **Summary:** The 2D Ising model is the minimal system that exhibits a phase transition ($T_c \approx 2.269 J/k_B$) and **spontaneous symmetry breaking** from simple, local interactions. Its exponentially large state space necessitates stochastic sampling.

#### Quiz Questions

**1. Which of the following phenomena is considered an **emergent property** of the Ising model?**

* **A.** The total number of spins on the lattice.
* **B.** The spontaneous symmetry breaking below $T_c$ (i.e., non-zero magnetization when $H=0$). (**Correct**)
* **C.** The value of the nearest-neighbor coupling $J$.
* **D.** The ability to write the Hamiltonian.

**2. The one-dimensional Ising model is analytically known to exhibit:**

* **A.** A sharp phase transition at a known $T_c$.
* **B.** A continuous phase transition with non-analytic behavior.
* **C.** **No phase transition** at any non-zero temperature ($T>0$). (**Correct**)
* **D.** Only antiferromagnetic ordering.



#### Interview-Style Question

**Question:** The Ising model is often called the "Hello, World!" of complex systems. Why do we study it today, given that it's an extreme simplification of a real magnet?

**Answer Strategy:** We study the Ising model primarily because it is the simplest system that rigorously demonstrates three crucial concepts:
1.  **Emergence:** Macroscopic collective order (magnetization) arising solely from simple local rules (nearest-neighbor interaction $J$).
2.  **Phase Transition/Criticality:** It is analytically solvable in 2D, providing an exact benchmark ($T_c$) to validate complex simulation methods like MCMC.
3.  **Universality:** The model's critical exponents describe an entire **universality class**, meaning its behavior near $T_c$ is shared by many real-world systems, regardless of their microscopic details.



### 2.2 The Hamiltonian and the Local Rule

> **Summary:** The energy (Hamiltonian) is defined by nearest-neighbor coupling $J$ and an external field $H$. For MCMC efficiency, the **change in energy $\Delta E$** for a single spin flip is calculated locally using only the spin and its four neighbors, making the operation $\mathcal{O}(1)$.

#### Quiz Questions

**1. For a **ferromagnetic** Ising model with coupling $J>0$, the Hamiltonian $E(\sigma)$ is minimized when neighboring spins are:**

* **A.** Anti-aligned ($\sigma_i \sigma_j = -1$).
* **B.** Uncorrelated.
* **C.** **Aligned** ($\sigma_i \sigma_j = +1$). (**Correct**)
* **D.** Coupled to a large external field $H$.

**2. The crucial computational advantage of using single-spin-flip Metropolis updates on the Ising model is that the calculation of $\Delta E$ is $\mathcal{O}(1)$. This is because:**

* **A.** $\Delta E$ is always zero.
* **B.** We recalculate the full energy $E(\sigma')$ and subtract $E(\sigma)$.
* **C.** Flipping a single spin only affects the energy contributions of that spin and its nearest neighbors. (**Correct**)
* **D.** The time step is very small.



#### Interview-Style Question

**Question:** In the context of the MCMC acceptance rule, $\alpha = \min(1, e^{-\beta \Delta E})$, describe the two scenarios for $\Delta E$ and what they tell us about the physics of the move.

**Answer Strategy:**
1.  **$\Delta E \le 0$ (Energy-Lowering or Neutral Move):** This move is **always accepted** ($\alpha=1$). This embodies the physical tendency of the system to seek the lowest energy state, quickly flowing down the energy landscape.
2.  **$\Delta E > 0$ (Energy-Increasing Move):** This move is **accepted with probability $e^{-\beta \Delta E}$**. This is the **thermal fluctuation** mechanism. At low temperature (large $\beta$), this probability is tiny, and the move is usually rejected. At high temperature (small $\beta$), the probability is near 1, allowing the system to easily overcome energy barriers, which is essential for ergodicity and exploring the state space.



### 2.3 The Computational Framework: Periodic Boundary Conditions

> **Summary:** We use **Periodic Boundary Conditions (PBCs)** to eliminate unphysical **surface effects** caused by edge spins having fewer neighbors than bulk spins. PBCs wrap the lattice onto a torus, ensuring every spin has the full coordination number (four in 2D) and preserving translational invariance.

#### Quiz Questions

**1. The primary purpose of using **Periodic Boundary Conditions** in an Ising simulation is to:**

* **A.** Prevent any energy-increasing spin flips.
* **B.** Ensure the Monte Carlo simulation runs in parallel.
* **C.** **Minimize finite-size and surface effects** by giving every spin the same number of neighbors. (**Correct**)
* **D.** Automatically calculate the correlation length.

**2. Implementing PBCs in code for a 2D lattice of size $N$ involves using which mathematical operation when calculating neighbor indices?**

* **A.** Multiplication.
* **B.** Division.
* **C.** The **modulo operator** (`%`). (**Correct**)
* **D.** The power function.



#### Interview-Style Question

**Question:** If you ran an Ising simulation with **open (free)** boundary conditions instead of periodic boundary conditions, how would this affect your measurement of the equilibrium magnetization $\langle |M| \rangle$ for a small $10 \times 10$ lattice?

**Answer Strategy:** On a small lattice, the surface-to-volume ratio is high.
* With open boundaries, spins at the edges have fewer stabilizing neighbors, making them **more susceptible to flipping**.
* This increased thermal fluctuation near the boundaries tends to **disorder** the system more easily than the bulk.
* The overall measured magnetization $\langle |M| \rangle$ for the entire lattice would therefore be **lower** than the true thermodynamic value, and the critical transition would appear less sharp or more **rounded**.



### 2.4 Implementing the Metropolis Algorithm on the Lattice

> **Summary:** A **Monte Carlo Sweep (MCS)** consists of $N^2$ attempted single-spin updates, ensuring every spin is considered once on average. We must manage **equilibration** and choose between the simpler **Metropolis** algorithm and the potentially faster **Heat-bath (Glauber)** algorithm.

#### Quiz Questions

**1. In an $N \times N$ Ising simulation, a single **Monte Carlo Sweep (MCS)** is defined as:**

* **A.** Running the simulation until the energy stabilizes.
* **B.** A single attempt to flip a spin.
* **C.** **$N^2$ successive single-spin update attempts**. (**Correct**)
* **D.** The time required to compute the autocorrelation function.

**2. Compared to the Metropolis method, the **Heat-bath (Glauber) dynamics** update is distinct because:**

* **A.** It only allows energy-lowering moves.
* **B.** It always results in a lower autocorrelation time.
* **C.** The new spin state is **directly sampled** from its local conditional probability, rather than accepting/rejecting a proposed flip. (**Correct**)
* **D.** It violates the detailed balance condition.



#### Interview-Style Question

**Question:** Your simulation is running very slowly, with a core loop only achieving 1 million spin flips per second. Suggest two simple coding/optimization steps that can be taken to significantly increase the performance of the local update loop.

**Answer Strategy:** The Metropolis core loop is dominated by computing $\Delta E$ and checking the exponential. Two primary optimization steps are:
1.  **Precomputing the Boltzmann Factors:** The local energy change $\Delta E$ can only take a small, finite set of values (e.g., $\pm 8J, \pm 4J, 0$). Pre-calculate and **tabulate** the acceptance probabilities $\exp(-\beta \Delta E)$ for all possible $\Delta E$ values to avoid costly $\exp()$ calls inside the inner loop.
2.  **Using Integer Spin Representation:** Representing spins as integer $\pm 1$ values (instead of floats or a boolean) allows the nearest-neighbor interactions to be computed using faster **integer arithmetic**.



### 2.5 Analysis I: Macroscopic Observables

> **Summary:** Macroscopic state is described by **observables** calculated from the microstates. The **Order Parameter** is the magnetization per spin $\langle |M| \rangle$, which is non-zero below $T_c$. The **Susceptibility $\chi$** and **Specific Heat $C_v$** are found by measuring the fluctuations of $M$ and $E$, respectively.

#### Quiz Questions

**1. Below $T_c$, the ferromagnetic Ising model exhibits **spontaneous ordering**. To detect this in a finite simulation, we must measure the:**

* **A.** Energy per spin $\langle e \rangle$.
* **B.** **Absolute magnetization per spin $\langle |M| \rangle$**. (**Correct**)
* **C.** Spin-spin correlation function $C(\mathbf{r})$.
* **D.** Variance of the kinetic energy.

**2. The **Specific Heat** $C_v$ is calculated in MCMC simulations by measuring the fluctuations of which microscopic quantity?**

* **A.** Magnetization $M$.
* **B.** **Total Energy $E$**. (**Correct**)
* **C.** Temperature $T$.
* **D.** Spin density $\rho$.



#### Interview-Style Question

**Question:** Explain the physical significance of the **magnetic susceptibility $\chi$** diverging (or peaking sharply in a finite system) exactly at the critical temperature $T_c$.

**Answer Strategy:** Susceptibility $\chi$ measures the system's **response to an external magnetic field $H$**. The formula $\chi \propto (\langle M^2 \rangle - \langle M \rangle^2)$ shows it is proportional to the **variance (fluctuations) of the magnetization**. At $T_c$:
* The correlation length diverges, meaning spins are highly correlated over the entire lattice.
* This critical state means the system is extremely sensitive to external perturbations. A tiny change in the external field $H$ can induce a massive change in the total magnetization $M$.
* The peak in $\chi$ is thus the computational signature of this highly unstable, highly correlated state right at the critical point.



### 2.6 Analysis II: Equilibration and Autocorrelation

> **Summary:** MCMC data is correlated along the Markov chain. **Thermalization (burn-in)** is the initial phase where the system relaxes to equilibrium, and all measurements must be discarded. **Autocorrelation** means successive measurements are not independent; this inflates error bars. This is corrected by estimating the **integrated autocorrelation time $\tau_{\text{int}}$** and using **data binning**.

#### Quiz Questions

**1. The practice of running a simulation for a period and **discarding** all initial data is known as:**

* **A.** Subsampling.
* **B.** **Thermalization** (or burn-in). (**Correct**)
* **C.** Critical slowing down.
* **D.** Finite-size scaling.

**2. If an observable's **integrated autocorrelation time $\tau_{\text{int}}$** is 50 sweeps, what does this tell us about the sampled data?**

* **A.** The system is out of equilibrium.
* **B.** The simulation must run for at least 50 sweeps.
* **C.** You need $\sim 50$ sweeps between measurements to get statistically independent samples. (**Correct**)
* **D.** The error is proportional to $\sqrt{50}$.



#### Interview-Style Question

**Question:** The simulation is run at $T=1.0$ (low temperature) and $T=4.0$ (high temperature). Explain why the **thermalization phase** might be significantly longer at the low temperature ($T=1.0$).

**Answer Strategy:** At low temperature, the system is dominated by ferromagnetic coupling $J$. This creates **large, stable domains** of aligned spins.
* To equilibrate from a random initial state, the system must form these large domains.
* The Metropolis acceptance probability for flipping a spin **inside** a large, stable domain is very low (high $\Delta E$), making the process of domain formation slow and difficult.
* At high temperature ($T=4.0$), thermal fluctuations are so large that domain boundaries are unstable, and the system rapidly disorders, leading to a much shorter thermalization time.



### 2.7 Core Application: Locating the Phase Transition

> **Summary:** The critical temperature $T_c$ is located by observing the sharp drop in $\langle |M| \rangle$ and the peaks in $\chi$ and $C_v$. The most accurate method for estimating the thermodynamic $T_c$ is by finding the intersection of the **Binder Cumulant $U_L(T)$** curves for different system sizes $L$.

#### Quiz Questions

**1. Which plot should be used to systematically locate the critical temperature $T_c$ in the thermodynamic limit ($L \to \infty$)?**

* **A.** $\langle E \rangle$ vs. $T$.
* **B.** $\langle |M| \rangle$ vs. $T$.
* **C.** **The Binder Cumulant $U_L(T)$** vs. $T$ for several $L$, and finding their intersection point. (**Correct**)
* **D.** The autocorrelation time $\tau_{\text{int}}$ vs. $T$.

**2. At the critical temperature $T_c$, the correlation length $\xi(T)$ is theoretically expected to:**

* **A.** Be zero.
* **B.** Diverge (become infinite). (**Correct**)
* **C.** Be equal to the lattice size $N$.
* **D.** Be exactly $2.269$.



#### Interview-Style Question

**Question:** You run a simulation and find that the peak in $\chi(T)$ is located at $T_{\text{peak}}=2.30$ for $L=32$ and $T_{\text{peak}}=2.28$ for $L=64$. Explain why the peak location shifts with $L$ and what the thermodynamic $T_c$ likely is.

**Answer Strategy:** The shift is a **finite-size effect**. In a finite system, the true critical behavior is **rounded off** because the correlation length $\xi$ cannot exceed the lattice size $L$. The peak in $\chi$ occurs when $\xi \approx L$. As the lattice size $L$ increases, the peak sharpens and shifts closer to the true, analytic critical temperature $T_c$. Since the exact $T_c \approx 2.269 J/k_B$, the shift suggests the thermodynamic critical point is at or very near the analytic value, and the $L=64$ result is a better approximation than $L=32$.



## 💡 Hands-On Simulation Projects (Chapter Conclusion)

These projects are designed to implement and test the core concepts of the Ising model, from the local update rule to the detection of the phase transition.

### Project 1: The Local Metropolis Update Rule (The Engine)

* **Goal:** Implement the core local update and $\mathcal{O}(1)$ energy calculation.
* **Setup:** Initialize a small $10 \times 10$ lattice with $\sigma_i = +1$ (all up). Use $J=1$ and $H=0$.
* **Steps:**
    1.  Write a function `calculate_delta_E(lattice, i, j)` that computes $\Delta E$ for flipping spin $(i, j)$ by looking only at its four neighbors and applying PBCs.
    2.  Write the Metropolis function `attempt_flip(lattice, i, j, beta)` that uses the calculated $\Delta E$ and the acceptance ratio.
    3.  Run a few thousand Monte Carlo sweeps at a very high $\beta$ (low $T$, e.g., $\beta=1.0$) and a very low $\beta$ (high $T$, e.g., $\beta=0.1$).
* ***Goal***: Confirm that the high-$\beta$ run mostly remains $\sigma_i=+1$ (low energy), while the low-$\beta$ run quickly becomes randomized (disordered).

### Project 2: Simulating the Magnetization Curve $\langle |M| \rangle(T)$

* **Goal:** Generate the classic S-shaped magnetization curve that reveals the phase transition.
* **Setup:** Use a fixed lattice size ($L=32$ or $L=64$), $J=1$, $H=0$.
* **Steps:**
    1.  Choose a temperature range $T \in [1.0, 4.0]$ (or $\beta \in [0.25, 1.0]$) with $\Delta T = 0.1$.
    2.  For each $T$:
        * Run a thermalization phase (e.g., $1000$ MCS) and discard data.
        * Run a measurement phase (e.g., $5000$ MCS) and record $|M|$ at each step.
        * Calculate the ensemble average $\langle |M| \rangle$ and $\langle M^2 \rangle$.
    3.  Plot $\langle |M| \rangle$ vs. $T$.
* ***Goal***: Visually identify the sharp drop near $T \approx 2.269$ and show the saturation at $\langle |M| \rangle \approx 1$ at low $T$ and $\langle |M| \rangle \approx 0$ at high $T$.

### Project 3: Visualizing Thermalization and Autocorrelation

* **Goal:** Quantify statistical error and justify the discarding of the burn-in phase.
* **Setup:** Run a simulation at the critical temperature $T_c \approx 2.269$ for $20,000$ sweeps.
* **Steps:**
    1.  Plot the raw time series of the Energy $E(t)$ for the full $20,000$ sweeps. Visually identify the burn-in period.
    2.  Compute and plot the **Autocorrelation Function $C_E(\tau)$** of the energy measurements (after removing burn-in).
    3.  Use the raw variance and the integrated autocorrelation time $\tau_{\text{int}}$ (summing $C_E(\tau)$) to calculate the statistically correct standard error of the mean $\text{Error} \propto \sqrt{2 \tau_{\text{int}} / N_{\text{meas}}}$.
* ***Goal***: Demonstrate the exponential decay of $C_E(\tau)$ and calculate the necessary spacing between measurements required to achieve reliable statistics.

### Project 4: Finding the Critical Exponent $\gamma$ (Advanced)

* **Goal:** Use the Susceptibility peak and finite-size scaling to confirm the critical behavior of the 2D Ising model.
* **Setup:** Run simulations for three different lattice sizes: $L=32$, $L=64$, and $L=128$.
* **Steps:**
    1.  For each $L$, sweep a fine temperature range around $T_c$ (e.g., $T \in [2.0, 2.5]$ with $\Delta T = 0.01$).
    2.  Measure $\langle M^2 \rangle$ and $\langle |M| \rangle$ to calculate the susceptibility $\chi_L(T) \propto (\langle M^2 \rangle - \langle |M| \rangle^2)$.
    3.  Plot $\chi_L(T)$ vs. $T$ for all three $L$ values. Observe the peaks getting taller and sharper as $L$ increases.
    4.  The scaling hypothesis states that $\chi_{\text{peak}} \propto L^{\gamma/\nu}$. Plot $\log(\chi_{\text{peak}})$ vs. $\log(L)$ and use linear regression to determine the slope $y = \gamma / \nu$.
* ***Goal***: Estimate the ratio $\gamma/\nu$ (which is analytically $1.75$ for the 2D Ising model) and show that MCMC can extract universal critical exponents.
