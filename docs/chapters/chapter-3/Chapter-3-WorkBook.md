## Chapter 3: Physics II: Lattice Gauge Theory (Workbook)

The goal of this chapter is to apply the MCMC engine (from Chapter 1) to a non-perturbative **Quantum Field Theory (QFT)**, showing how a discrete lattice simulation can solve problems unreachable by continuous mathematical expansion.

| Section | Topic Summary |
| :--- | :--- |
| **3.1** | The Non-Perturbative Universe: Why We Need LGT |
| **3.2** | The Field on the Lattice: Sites, Links, and Plaquettes |
| **3.3** | The Wilson Action (The LGT "Hamiltonian") |
| **3.4** | The Simulation: Metropolis Updates on $SU(N)$ Matrices |
| **3.5** | Core Application: Measuring Confinement with the Wilson Loop |

---

### 3.1 The Non-Perturbative Universe

> **Summary:** We need LGT because the strong nuclear force (QCD) is "strongly coupled" at low energies, making traditional **perturbation theory** fail. The **Wick rotation** transforms the path integral into a statistical mechanics partition function, enabling MCMC simulation.

#### Quiz Questions

**1. The primary reason traditional perturbation theory fails in Quantum Chromodynamics (QCD) at low energies is:**

* **A.** The coupling constant $g$ is too small to be used for expansion.
* **B.** The coupling constant $g$ is too large (strong coupling), and the theory exhibits non-perturbative phenomena like confinement. (**Correct**)
* **C.** The QCD action is non-gauge-invariant.
* **D.** All gauge theories must be simulated on a lattice.

**2. The purpose of the **Wick rotation** in Lattice Gauge Theory (LGT) is to:**

* **A.** Replace the spatial dimensions with temporal dimensions.
* **B.** Convert the oscillatory factor $\exp(iS/\hbar)$ in the path integral into a real, non-oscillatory Boltzmann-like weight $\exp(-S_E/\hbar)$. (**Correct**)
* **C.** Force all degrees of freedom to be matrices instead of scalar fields.
* **D.** Ensure the lattice spacing $a$ goes to zero.



#### Interview-Style Question

**Question:** Conceptually, LGT is often described as a technique for studying **emergent phenomena**. Briefly explain what "emerges" in LGT (in the context of the strong force) that is not visible in the elementary equations.

**Answer Strategy:** The fundamental equations describe quarks and gluons. The emergent phenomenon is **confinement**. The emergence of a linearly rising potential $V(R) \sim \sigma R$ (the Area Law) and the physical observation of colour-neutral **hadrons** (like protons), rather than free quarks, are not directly predicted by the simple perturbative equations. They arise non-perturbatively from the collective, strong-coupling dynamics of the gauge fields sampled by MCMC.



### 3.2 The Field on the Lattice: Sites, Links, and Plaquettes

> **Summary:** The 4D spacetime is discretized. Gauge fields (gluons) are placed on the **links** as $SU(N)$ matrices. Matter fields (quarks) live on the **sites**. The fundamental gauge-invariant building block is the **plaquette** (a $1 \times 1$ closed loop of links).

#### Quiz Questions

**1. In Lattice Gauge Theory, the fundamental dynamical variables that we update using MCMC are the:**

* **A.** Scalar fields living on the sites.
* **B.** $SU(N)$ matrices living on the links, representing parallel transport. (**Correct**)
* **C.** Energy values assigned to each plaquette.
* **D.** The field strength tensor $F_{\mu\nu}$ at each site.

**2. Which object is the *smallest* and most fundamental structure on the lattice that is guaranteed to be **gauge-invariant** (up to conjugation)?**

* **A.** A single site.
* **B.** A single link $U_\mu(n)$.
* **C.** The trace of the plaquette product $\mathrm{Tr}\,U_{\mu\nu}(n)$. (**Correct**)
* **D.** The product of two adjacent links.



#### Interview-Style Question

**Question:** Why is it necessary to put the gauge field variables (gluons) on the **links** between sites rather than directly on the **sites** like a traditional scalar field?

**Answer Strategy:** This is required to maintain **local gauge invariance**. The link variable $U_\mu(n)$ represents the operation of **parallel transport**—a fundamental concept in gauge theory. When a gauge transformation $\Omega(n)$ is applied at a site, the link variable must be rotated by $\Omega(n)$ at its tail and inversely rotated by $\Omega^\dagger(n+\hat{\mu})$ at its head. Placing the field on the link ensures this property is preserved exactly on the lattice, coupling the degrees of freedom consistently between neighboring sites.



### 3.3 The Wilson Action (The LGT “Hamiltonian”)

> **Summary:** The continuous Yang-Mills action is approximated by the **Wilson Action** $S_W$, which is a sum over the traces of all plaquettes. The coupling parameter $\beta_{\text{gauge}} = 2N/g_0^2$ acts like the **inverse temperature** $\beta$ from the Ising model, controlling the transition between the strong-coupling (confining) and weak-coupling (deconfining) regimes.

#### Quiz Questions

**1. The bare lattice coupling $\beta_{\text{gauge}}$ in the Wilson Action is directly analogous to which physical parameter in the Ising model?**

* **A.** The coupling constant $J$.
* **B.** The external magnetic field $H$.
* **C.** The inverse temperature $\beta = 1/k_B T$. (**Correct**)
* **D.** The system size $N \times N$.

**2. The Wilson Action $S_W$ is designed to be minimized when the plaquette matrix $U_{\mu\nu}(n)$ is close to the:**

* **A.** Zero matrix.
* **B.** Identity matrix $I$. (**Correct**)
* **C.** Negative of the identity matrix $-I$.
* **D.** The Pauli matrices $\sigma^a$.



#### Interview-Style Question

**Question:** In the context of LGT, how does the physicist ensure that their simulated lattice results, obtained at a finite lattice spacing $a$, are physically meaningful for the continuous world? What parameter must be controlled?

**Answer Strategy:** To recover physical meaning, the simulation must be performed in the **continuum limit**, where the lattice spacing $a$ approaches zero. This is achieved by tuning the bare lattice coupling $\beta_{\text{gauge}}$. The continuum limit corresponds to taking $\beta_{\text{gauge}} \to \infty$ (weak coupling), where the lattice artifacts are minimized, and the physical results become independent of the specific lattice parameter $a$. This scaling is essential for all LGT calculations.



### 3.4 The Simulation: Metropolis Updates on $SU(N)$ Matrices

> **Summary:** Metropolis updates in LGT use a **local update** strategy. The change in action $\Delta S$ depends only on the local **staple** $\Sigma_\mu(n)$—the products of the three neighboring links in the plaquettes connected to the updated link. This $O(1)$ calculation is crucial for efficiency. Updates must use a random matrix $\delta U$ that preserves the $SU(N)$ group manifold.

#### Quiz Questions

**1. When updating a single link matrix $U_\mu(n)$, the change in the action $\Delta S_W$ depends only on the sum of the products of the *three* neighboring links in the surrounding plaquettes. This sum is formally called the:**

* **A.** Wilson loop.
* **B.** Staple $\Sigma_\mu(n)$. (**Correct**)
* **C.** Propagator.
* **D.** Adjacency matrix.

**2. An efficient Metropolis update in LGT requires proposing a new link matrix $U'$ that is:**

* **A.** An arbitrary random matrix.
* **B.** Drawn from a uniform distribution over the entire $SU(N)$ group manifold.
* **C.** Close to the current matrix $U$ to maintain a reasonable acceptance rate. (**Correct**)
* **D.** Always accepted (acceptance probability $\alpha=1$).



#### Interview-Style Question

**Question:** Contrast the **Metropolis** update with the **Heat-Bath** update in LGT, and explain why both are often preferred over a third method, **Over-relaxation**, for generating statistically independent configurations.

**Answer Strategy:**
* **Metropolis** and **Heat-Bath** both satisfy **detailed balance** and are used to generate samples whose distribution converges to the target $\exp(-S_W)$. Heat-bath is generally more efficient for $SU(2)$ as it generates the new configuration *exactly* from the conditional probability.
* **Over-relaxation** is a *deterministic* update that preserves the action (microcanonical). While it reduces autocorrelation, it does **not** explore new energy surfaces. Therefore, it must be *combined* with Metropolis or Heat-Bath steps to ensure **ergodicity** (the ability to reach all relevant configurations) and maintain the correct canonical distribution.



### 3.5 Core Application: Measuring Confinement with the Wilson Loop

> **Summary:** The **Wilson Loop $W(R, T)$** (a product of link matrices around an $R \times T$ rectangle) is the key observable to test confinement. Confinement implies a linearly rising potential $V(R) \sim \sigma R$, which computationally translates to the **Area Law**: $\langle W(R, T) \rangle \sim e^{-\sigma R T}$, where the decay rate is proportional to the loop's area $R \times T$.

#### Quiz Questions

**1. The primary physical phenomenon that the **Area Law** behavior of the Wilson Loop demonstrates in pure QCD is:**

* **A.** Asymptotic freedom at high energy.
* **B.** The breaking of chiral symmetry.
* **C.** Quark **confinement**. (**Correct**)
* **D.** The existence of glueballs.

**2. The expectation value $\langle W(R,T) \rangle$ is computationally related to the static quark-antiquark potential $V(R)$ by:**

* **A.** $\langle W(R,T) \rangle \propto V(R)/T$.
* **B.** $\langle W(R,T) \rangle \propto e^{-V(R) T}$. (**Correct**)
* **C.** $\langle W(R,T) \rangle \propto V(R) / R$.
* **D.** $\langle W(R,T) \rangle \propto R \times T$.



#### Interview-Style Question

**Question:** The Wilson Loop observable suffers from a poor signal-to-noise ratio at large areas. Explain why this happens and propose one computational technique to mitigate this issue.

**Answer Strategy:** The signal $\langle W(R, T) \rangle$ decays **exponentially** with the loop area $R \times T$. However, the noise (the variance of the product of link matrices) does not decrease at the same rate. This means for large loops, the signal is drowned out by noise.
One mitigation technique is **link smearing** (e.g., APE smearing), which replaces a link with a weighted average of itself and its surrounding links (the staple). This process smooths out short-distance, high-frequency fluctuations (noise) without changing the physics, thus enhancing the signal's visibility and making large loops tractable.



### Hands-On Project: The Area Law as an Exponential Fit 🧪

This project simulates the key analytical step required to prove confinement, using conceptual data.

* **Project Goal:** Write a numerical routine to extract the String Tension ($\sigma$) by fitting the Wilson loop data to the Area Law, demonstrating the exponential decay.
* **Setup:**
    1.  Define a target string tension, e.g., $\sigma_{\text{target}} = 0.05$.
    2.  Generate a set of $R \times T$ areas, $A_k = R_k \times T_k$, for $k=1, \dots, 20$.
    3.  Generate "simulated" Wilson Loop expectation values $W_k$ that follow the Area Law, but include small random Gaussian noise $\eta_k$:
        $$\mathbf{W}_k = \exp(-\sigma_{\text{target}} \cdot A_k) + \eta_k$$
* **Computational Strategy:**
    1.  Transform the equation into a linear form: $\ln(\mathbf{W}_k) \approx - \sigma_{\text{target}} \cdot A_k$.
    2.  Use a linear least-squares fitting routine (e.g., `np.polyfit(A_k, np.log(W_k), 1)`) to fit the data.
    3.  Extract the slope of the linear fit, which is the estimated string tension $\sigma_{\text{fit}}$.
* **Visualization:** Plot the raw data ($\ln(\mathbf{W}_k)$ vs. $A_k$) and the best-fit line.
* ***Expected Result:*** The plot will show a clear linear decay, and the fitted slope $\sigma_{\text{fit}}$ will be close to the target value $\sigma_{\text{target}}$, confirming the computational proof of the Area Law.


