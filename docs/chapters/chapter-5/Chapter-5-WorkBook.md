## Chapter 5: Biology I: Stochastic Systems Biology (Workbook)

The goal of this chapter is to shift from continuous stochastic models (like GBM) to **discrete stochastic models**, introducing the **Gillespie Algorithm (SSA)** as the mathematically exact method for simulating biochemical reactions dominated by low copy number noise.

| Section | Topic Summary |
| :--- | :--- |
| **5.1** | Chapter Opener: Noise is the Whole Story |
| **5.2** | The Gillespie Algorithm (SSA) |
| **5.3** | Simulation: Modeling Simple Gene Expression |
| **5.4** | Application: Observing Transcriptional Bursting |
| **5.5** | Chapter Summary & Bridge to Chapter 6 |



### 5.1 Noise is the Whole Story

> **Summary:** Gene expression is a stochastic process where discrete, random events lead to significant **cell-to-cell variability** (noise). In low copy number regimes, this noise is significant, necessitating stochastic models over traditional deterministic ODEs.

#### Section Detail

Noise is classified into **intrinsic noise** (randomness within a single cell, leading to uncorrelated fluctuations in identical genes) and **extrinsic noise** (cell-to-cell variability in global factors, leading to correlated fluctuations). Noise is quantified using the **Fano Factor ($\eta$)**, which is the ratio of variance to mean ($\eta = \mathrm{Var}(X) / \langle X\rangle$). Super-Poissonian noise ($\eta > 1$) is often the signature of **transcriptional bursting**.

#### Quiz Questions

**1. Which dimensionless measure is used to quantify the strength of gene expression noise by comparing the variance to the mean?**

* **A.** The Poisson Ratio.
* **B.** The Fano Factor ($\eta = \mathrm{Var}(X)/\langle X\rangle$). (**Correct**)
* **C.** The Autocorrelation Function.
* **D.** The $\mathcal{O}(1)$ complexity factor.

**2. **Intrinsic Noise** in gene expression arises primarily from:**

* **A.** Cell-to-cell differences in cell size or cell cycle stage (extrinsic factors).
* **B.** **The random timing and nature of biochemical reactions** within a single cell. (**Correct**)
* **C.** Large copy numbers of mRNA molecules.
* **D.** Transcription and translation occurring at a perfectly steady rate.



#### Interview-Style Question

**Question:** An experiment measures the concentration of a protein and finds its Fano Factor is $\eta = 5.0$. Explain what this result suggests about the underlying molecular production process compared to a simple, continuous-rate birth-and-death process.

**Answer Strategy:** A simple, continuous birth-and-death process (like a Poisson process) should yield a Fano Factor of $\eta=1$. Finding $\eta = 5.0$ (which is $\eta > 1$, or **super-Poissonian noise**) strongly suggests that the production is **bursty**. This means that instead of producing molecules at a steady, continuous rate, the gene produces them in large, episodic bursts separated by long periods of inactivity, leading to large fluctuations in copy number.



### 5.2 The Gillespie Algorithm (SSA)

> **Summary:** The **Gillespie Algorithm (SSA)** provides an exact procedure for simulating the discrete trajectories consistent with the **Chemical Master Equation (CME)**. It uses two random numbers to determine the **time $\tau$ to the next reaction** and **which reaction $j$ occurs**, both governed by the **propensity functions $a_j(\mathbf{x})$**.

#### Section Detail

The SSA is a continuous-time Monte Carlo method. Its exactness relies on the Poisson process assumption for reaction events.
1.  **Waiting Time ($\tau$):** The time until the next reaction is exponentially distributed with parameter $a_0(\mathbf{x}) = \sum_j a_j(\mathbf{x})$ (the total propensity): $$\tau = \frac{1}{a_0} \ln\left(\frac{1}{r_1}\right)$$
2.  **Reaction Choice ($j$):** The probability of reaction $j$ occurring is proportional to its propensity: $\mathbb{P}(j) = a_j(\mathbf{x}) / a_0(\mathbf{x})$.

This procedure generates a statistically exact trajectory of molecule counts over time.

#### Quiz Questions

**1. In the Gillespie Algorithm, the waiting time ($\tau$) until the next reaction is drawn from which probability distribution?**

* **A.** A standard Normal distribution.
* **B.** A discrete binomial distribution.
* **C.** An **Exponential distribution**. (**Correct**)
* **D.** A uniform distribution.

**2. Which mathematical equation governs the time evolution of the probability distribution over the vector of molecule counts in a well-mixed stochastic system?**

* **A.** The Black-Scholes-Merton Equation.
* **B.** The Schrödinger Equation.
* **C.** The **Chemical Master Equation (CME)**. (**Correct**)
* **D.** The Newton-Raphson formula.



#### Interview-Style Question

**Question:** The computational cost of the basic SSA scales linearly with the number of reaction channels, $\mathcal{O}(M)$. Why does the algorithm have to recompute and re-sum all $M$ reaction propensities at every single time step?

**Answer Strategy:** The SSA must recompute all propensities because the fundamental input for the next step—the **total propensity $a_0(\mathbf{x})$**—changes every time a reaction occurs. Since $a_0(\mathbf{x})$ defines the rate of the Poisson process, it changes the mean waiting time ($\tau = 1/a_0$). If $a_0$ were held constant, the simulation would violate the memoryless property of the Markov chain. More efficient variants exist (e.g., Next Reaction Method) that use priority queues to reduce this overhead to $\mathcal{O}(\log M)$.



### 5.3 Simulation: Modeling Simple Gene Expression

> **Summary:** The minimal gene expression model (the **telegraph model**) comprises gene switching (ON/OFF), transcription ($m$ production), and molecule degradation ($\gamma_m m$). Propensity functions are written to reflect current state, such as $a_{\text{transcription}} = k_m g$, where $g \in \{0, 1\}$ is the gene state.

#### Section Detail

The telegraph model is simplified to capture the core of bursting.
* Gene activation and inactivation rates ($k_{\text{on}}, k_{\text{off}}$) determine the burst frequency.
* Transcription ($k_m$) and mRNA degradation ($\gamma_m$) rates determine the burst size, with mean burst size being $k_m / \gamma_m$.
The simulation generates a time series of molecule counts, $\mathbf{x}(t)$, which, when run many times, yields the steady-state distribution.

#### Quiz Questions

**1. In the simple gene expression telegraph model, the propensity function for mRNA decay ($m \xrightarrow{\gamma_m} \emptyset$) is $a_4 = \gamma_m m$. This form signifies that:**

* **A.** Decay occurs only when the gene is ON.
* **B.** The likelihood of decay is proportional to the number of existing mRNA molecules ($m$). (**Correct**)
* **C.** Decay is a second-order (bimolecular) reaction.
* **D.** The rate of decay is constant and independent of the state.

**2. The expected **mean burst size** (average number of mRNAs produced in one burst) in the telegraph model is determined by the ratio of which two rates?**

* **A.** $k_{\text{on}} / k_{\text{off}}$.
* **B.** $k_{\text{off}} / k_m$.
* **C.** **Transcription rate ($k_m$) / mRNA degradation rate ($\gamma_m$)**. (**Correct**)
* **D.** $k_m / k_{\text{off}}$.



#### Interview-Style Question

**Question:** An engineer wants to use the telegraph model to minimize protein noise. They decide to reduce the transcription rate $k_m$ by a factor of 10 and increase the translation rate $k_p$ by a factor of 10. Analyze how these changes would affect the noise characteristics ($\eta$) and the mean protein level.

**Answer Strategy:**
1.  **Mean Protein Level ($\langle p \rangle$):** The steady-state mean is proportional to the product of production rates ($k_m \cdot k_p$) and inversely proportional to degradation rates. Since $k_m$ decreases by 10 and $k_p$ increases by 10, the **mean protein level remains roughly unchanged** (assuming other parameters are constant).
2.  **Noise ($\eta$):** Noise is largely driven by burst size $k_m / \gamma_m$. Reducing $k_m$ by a factor of 10 means the mean **burst size is reduced by a factor of 10**. Smaller bursts smooth out fluctuations, meaning the system transitions to a less-bursty, less-noisy regime. This change should **reduce the overall Fano Factor** ($\eta$), making the expression less super-Poissonian.


### 5.4 Application: Observing Transcriptional Bursting

> **Summary:** SSA simulations generate trajectories that exhibit **transcriptional bursting**—episodic mRNA production—which explains the super-Poissonian noise ($\eta > 1$) observed experimentally. Simulations provide the full **distribution** of molecule counts, contrasting with deterministic ODEs that only yield the mean.

#### Section Detail

Bursting occurs because the gene stochastically switches between the ON and OFF promoter states. Analyzing the SSA trajectories allows us to estimate the burst frequency and size directly. The importance of the SSA lies in its ability to capture **phenotypic heterogeneity**—the cell-to-cell variability in molecule counts—which is crucial for understanding biological decisions but invisible to mean-field ODE models.

#### Quiz Questions

**1. When simulating the mRNA copy number over time using the Gillespie algorithm for the telegraph model, the typical trajectory visually resembles a:**

* **A.** Smooth, exponentially decaying curve.
* **B.** A deterministic, fixed-point equilibrium line.
* **C.** A **"spike train" or episodic bursts** of production followed by exponential decay. (**Correct**)
* **D.** A continuous log-normal random walk.

**2. A key failure of deterministic ODE models in the context of gene expression is that they cannot capture:**

* **A.** The average expression level.
* **B.** The fact that expression levels can be positive.
* **C.** The **probability distribution of expression levels (variance)** across a population of cells. (**Correct**)
* **D.** The steady-state mean concentration.



#### Interview-Style Question

**Question:** If you were simulating a batch of 100 identical cells using the SSA, and you notice that after 10,000 seconds, the fluctuations in all 100 cell trajectories look very similar (correlated), what source of noise would you investigate first?

**Answer Strategy:** Correlated fluctuations across a population of genetically identical cells suggest a high contribution of **Extrinsic Noise**. Extrinsic noise is due to cell-to-cell variability in global factors that affect *all* genes similarly, such as ribosome or polymerase copy numbers. The simulation would need to be modified by sampling these global parameters from a distribution *between* the 100 cell runs, rather than assuming they are identical.



## 💡 Hands-On Simulation Projects (Chapter Conclusion) 🧪

These projects are designed to implement the core Gillespie Algorithm and demonstrate its ability to capture molecular noise and bursting.

### Project 1: Implementing the Gillespie Core (The Engine)

* **Goal:** Implement the fundamental SSA step to generate the waiting time and select the reaction.
* **Setup:** Define a minimal system with two reactions: **A $\xrightarrow{k_1}$ B** and **B $\xrightarrow{k_2}$ $\emptyset$**, with initial counts $N_A=10, N_B=0$ and rates $k_1=1.0, k_2=0.1$.
* **Steps:**
    1.  Define the propensities $a_1 = k_1 N_A$ and $a_2 = k_2 N_B$.
    2.  Write a function `gillespie_step(NA, NB)` that:
        * Calculates $a_0 = a_1 + a_2$.
        * Draws $r_1$ and $r_2$.
        * Calculates $\tau$ and selects the reaction $j$ using the cumulative sum rule.
        * Returns $(\tau, j)$.
    3.  Run a loop for a few hundred steps, printing the time and $N_A, N_B$ after each update.
* ***Goal***: Verify that the simulation time advances stochastically and that reaction $j=1$ is 10 times more likely to occur than $j=2$ initially.

### Project 2: Simulating and Visualizing Transcriptional Noise

* **Goal:** Implement the full telegraph model (mRNA only) and compare the stochastic trajectory to the deterministic mean.
* **Setup:** Use the parameters: $k_{\text{on}}=0.01, k_{\text{off}}=0.1, k_m=1.0, \gamma_m=0.05$. Initial state $g=0, m=0$. Final time $T_{\text{final}} = 1000$ seconds.
* **Steps:**
    1.  Implement the full SSA loop for the gene switching, transcription, and mRNA decay reactions (4 channels).
    2.  Run the SSA once and plot the time series of the mRNA count $m(t)$.
    3.  Calculate the deterministic steady-state mean mRNA count: $\langle m \rangle_{\text{det}} = \frac{k_{\text{on}}}{k_{\text{on}} + k_{\text{off}}} \cdot \frac{k_m}{\gamma_m}$.
    4.  Plot the deterministic mean as a horizontal line over the stochastic trajectory.
* ***Goal***: Show that the stochastic trajectory exhibits large, irregular bursts around the smoother deterministic mean, visually demonstrating the noise.

### Project 3: Quantifying Super-Poissonian Noise ($\eta$)

* **Goal:** Quantitatively demonstrate that the telegraph model produces super-Poissonian noise (Fano Factor $\eta > 1$).
* **Setup:** Use the bursty parameters from Project 2.
* **Steps:**
    1.  Run $M=500$ independent SSA trajectories until they all reach the steady state ($T=10,000$ s).
    2.  Record the final mRNA count $m_i$ for each of the $M$ trajectories.
    3.  Calculate the ensemble mean $\langle m \rangle = \frac{1}{M} \sum m_i$ and the ensemble variance $\mathrm{Var}(m) = \frac{1}{M-1} \sum (m_i - \langle m \rangle)^2$.
    4.  Compute the Fano Factor $\eta = \mathrm{Var}(m) / \langle m \rangle$.
* ***Goal***: Show that the calculated $\eta$ is significantly greater than $1.0$ (e.g., $\eta \approx 10-20$), confirming the noise is dominated by large-size transcriptional bursts.
