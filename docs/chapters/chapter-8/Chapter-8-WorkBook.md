## 📉 Chapter 8: Finance II: The Stochastic Calculus (SDEs) (Workbook)

The goal of this chapter is to introduce the necessary mathematical framework — **Itō Calculus** — to model asset prices as continuous, random processes and to numerically solve their **Stochastic Differential Equations (SDEs)**.

| Section | Topic Summary |
| :--- | :--- |
| **8.1** | Chapter Opener: Why Classical Calculus Fails |
| **8.2** | The Foundation: The Wiener Process $W_t$ and SDEs |
| **8.3** | The Breakthrough: Itō’s Lemma |
| **8.4** | The Solver: The Euler–Maruyama Method |
| **8.5** | Chapter Summary & Bridge to Chapter 9 |

***

### 8.1 Why Classical Calculus Fails

> **Summary:** Classical calculus fails for asset prices because the price path is **continuous but nowhere differentiable**. This chaotic motion, which approximates the **Wiener Process ($W_t$)**, violates the smoothness assumption needed for standard differentiation.

#### Section Detail

Classical differentiation requires that the rate of change is well-behaved and deterministic. The core issue with $W_t$ is the unexpected scale of its quadratic variation: the square of its infinitesimal increment, $(dW_t)^2$, is of order $dt$, not $(dt)^2$. This means key terms that are dropped in classical Taylor expansions must be kept in the new **Itō calculus**.

#### Quiz Questions

**1. The failure of classical calculus to model stock prices is primarily due to the fact that the price path is:**

* **A.** A deterministic sine wave.
* **B.** A discrete random walk.
* **C.** **Continuous but nowhere differentiable**. (**Correct**)
* **D.** A simple linear function of time.

**2. Which of the following is the key reason why standard calculus rules (like dropping second-order terms in Taylor expansions) break down for the Wiener Process?**

* **A.** The process has infinite drift.
* **B.** **The square of the differential, $(dW_t)^2$, is of order $dt$, not $(dt)^2$**. (**Correct**)
* **C.** The process is not a Markov chain.
* **D.** The time step $\Delta t$ is always too large.

---

#### Interview-Style Question

**Question:** In simple terms, explain the mathematical significance of the relationship: $\mathbb{E}[(dW_t)^2] = dt$.

**Answer Strategy:** This equation is the core difference between classical and stochastic calculus. It means the variance of the infinitesimal price shock ($\mathbb{E}[(dW_t)^2]$) grows **linearly with time $dt$**. In classical calculus, any term of order $dt$ in a Taylor expansion is assumed negligible compared to $dt$. By showing that the variance of the random component is *also* of order $dt$, we demonstrate that **the random component is not negligible** and must be explicitly retained, leading to the necessary modification of the chain rule.

---

***

### 8.2 The Foundation: The Wiener Process $W_t$ and SDEs

> **Summary:** The **Wiener Process ($W_t$)** is the continuous limit of the random walk. It is defined by its **Gaussian increments** ($\sim \mathcal{N}(0, \Delta t)$). It drives the general **Stochastic Differential Equation (SDE)**: $dS_t = \mu\,dt + \sigma\,dW_t$, which combines deterministic **drift ($\mu$)** and stochastic **diffusion ($\sigma$)**.

#### Section Detail

The SDE is the mathematical tool for describing continuous random dynamics. The most important example in finance is **Geometric Brownian Motion (GBM)**, $dS_t = \mu S_t\,dt + \sigma S_t\,dW_t$, which models prices as log-normal processes. The term $\mu\,dt$ governs the long-term trend, while the $\sigma\,dW_t$ term governs the short-term volatility (or random shock).

#### Quiz Questions

**1. The SDE $dS_t = \mu(S_t,t)\,dt + \sigma(S_t,t)\,dW_t$ contains two parts. The term $\mu(S_t,t)\,dt$ is known as the:**

* **A.** Volatility term.
* **B.** **Drift term**. (**Correct**)
* **C.** Diffusion term.
* **D.** Stochastic noise.

**2. The single most widely used SDE for modeling asset prices in financial mathematics is:**

* **A.** The Mean-Reverting Ornstein-Uhlenbeck process.
* **B.** **Geometric Brownian Motion (GBM)**. (**Correct**)
* **C.** The Pure Diffusion equation.
* **D.** The Ito-Correction SDE.

---

#### Interview-Style Question

**Question:** In the context of GBM, $dS_t = \mu S_t\,dt + \sigma S_t\,dW_t$, what is the physical meaning of the diffusion term being proportional to the current price ($ \propto S_t$)?

**Answer Strategy:** This proportionality ensures that the price remains positive and reflects a key financial reality: **volatility scales with price**. A $\$100$ stock has a much larger dollar movement (volatility) than a $\$1$ stock. By making the diffusion term proportional to $S_t$, we ensure that the *percentage change* in the price remains constant (log-normal property), rather than the absolute dollar change, providing a more realistic model for financial markets.

---

***

### 8.3 The Breakthrough: Itō’s Lemma

> **Summary:** **Itō’s Lemma** is the stochastic analog of the chain rule. It extends the classical chain rule by adding the necessary **Itō correction term** $\frac{1}{2}\sigma^2 \frac{\partial^2 f}{\partial S^2}\,dt$, which accounts for the accumulated drift caused purely by volatility.

#### Section Detail

The need for Itō’s Lemma comes directly from the rule $(dW_t)^2 = dt$, which forces the second-order Taylor term $\frac{1}{2}\frac{\partial^2 f}{\partial S^2}(dS_t)^2$ to be retained and simplified to the deterministic correction. This correction ensures that when a function of a stochastic variable is differentiated, the result is consistent with the underlying randomness. The exact solution for GBM, which involves the $\mu - \tfrac{1}{2}\sigma^2$ term, is a direct result of applying Itō's Lemma to $\ln S_t$.

$$
df = \left(\frac{\partial f}{\partial t} + \mu \frac{\partial f}{\partial S} + \frac{1}{2}\sigma^2 \frac{\partial^2 f}{\partial S^2}\right) dt + \sigma \frac{\partial f}{\partial S}\,dW_t
$$

#### Quiz Questions

**1. The key non-classical rule that underlies the derivation of Itō’s Lemma is that:**

* **A.** $dt^2 = 0$.
* **B.** $dW_t\,dt = 0$.
* **C.** **$(dW_t)^2 = dt$**. (**Correct**)
* **D.** $\mu = r$.

**2. The term $\frac{1}{2}\sigma^2 \frac{\partial^2 f}{\partial S^2}\,dt$ in Itō’s Lemma is known as the Itō correction. It is fundamentally a:**

* **A.** Random, stochastic term.
* **B.** **Deterministic drift adjustment**. (**Correct**)
* **C.** Second-order noise term.
* **D.** First-order velocity term.

---

#### Interview-Style Question

**Question:** Consider the SDE for GBM. If you apply the classical chain rule to $f(S_t) = S_t^2$ and the Itō Lemma to $f(S_t) = S_t^2$, the two results differ by a term proportional to $\sigma^2 dt$. Why is the term $\propto \sigma^2 dt$ always present in the Itō version and missing in the classical version?

**Answer Strategy:** The term $S_t^2$ is convex ($\frac{\partial^2 f}{\partial S^2} > 0$).
* The **classical rule** ignores the second-order term $\frac{1}{2}f_{SS}(dS_t)^2$.
* The **Itō rule** keeps this term, which simplifies to $\frac{1}{2}(2)\sigma^2 dt = \sigma^2 dt$.
The missing term represents the deterministic drift that the price gains *due to its own volatility* ($\sigma$). Since the path is always jiggling (volatility $\sigma>0$), the function is always growing slightly faster than predicted by the average trend, and Itō's Lemma correctly captures this gain.

---

***

### 8.4 The Solver: The Euler–Maruyama Method

> **Summary:** The **Euler–Maruyama (EM) method** is the simplest numerical scheme for solving SDEs. It discretizes the SDE by taking a step proportional to the drift ($\Delta t$) and adding a random shock proportional to the diffusion ($\sqrt{\Delta t} Z$). Its **weak convergence** is $O(\Delta t)$, which is often sufficient for calculating expected payoffs in finance.

#### Section Detail

The EM formula is $S_{t+\Delta t} \approx S_t + \mu(S_t,t)\Delta t + \sigma(S_t,t)\sqrt{\Delta t}Z_t$, where $Z_t \sim \mathcal{N}(0,1)$.

| Type | Definition | Order of Convergence |
| :--- | :--- | :--- |
| **Strong Convergence** | Pathwise accuracy (individual trajectory) | $O(\sqrt{\Delta t})$ |
| **Weak Convergence** | Accuracy of expected values (mean) | $O(\Delta t)$ |

For options, we primarily need weak convergence, making EM the preferred, stable method.

#### Quiz Questions

**1. The primary random component added at each step of the Euler–Maruyama simulation must be scaled by:**

* **A.** The square of the time step $(\Delta t)^2$.
* **B.** The total time $T$.
* **C.** **The square root of the time step $\sqrt{\Delta t}$**. (**Correct**)
* **D.** The initial price $S_0$.

**2. In financial modeling, we often use the Euler–Maruyama method for its **weak convergence** because it accurately estimates the:**

* **A.** Exact pathwise solution of a single trajectory.
* **B.** Required computational time.
* **C.** **Expected value (average payoff) of the SDE**. (**Correct**)
* **D.** Strong order of convergence.

---

#### Interview-Style Question

**Question:** The numerical simulation of an SDE is required for Monte Carlo option pricing. Why is the simplicity and stability of the first-order **Euler–Maruyama** method often favored in financial practice over a more accurate higher-order method like Milstein?

**Answer Strategy:** In finance, the goal is typically to find the expected payoff, which relies on **weak convergence** ($O(\Delta t)$). Since EM already achieves a weak order of $O(\Delta t)$, the complexity added by higher-order methods (which offer marginal gains in weak accuracy but are more complex to implement and debug) is generally not worth the effort. The simple structure of EM and its stability make it the most reliable and transparent choice for production systems.

---

***

### 8.5 Chapter Summary & Bridge to Chapter 9

> **Summary:** Stochastic calculus provides the tools to model and simulate randomness. The **Itō Correction** reveals that volatility adds a predictable deterministic component to the drift. This prepares the way for the **Black–Scholes–Merton (BSM) derivation**, where the random terms are perfectly **cancelled out** in a hedged portfolio, resulting in a deterministic **PDE** for the option price.

#### Section Detail

Chapter 8 provided the mathematical foundation for Part II's finance section. Itō’s Lemma provides the analytical tool, and the EM method provides the numerical integration tool. The key takeaway for Chapter 9 is the recognition that the random term in the derivative's SDE (the $dW_t$ term) is directly proportional to the random term in the underlying asset's SDE. This symmetry allows for its complete cancellation via dynamic hedging.

#### Quiz Questions

**1. The philosophical leap achieved by the Black–Scholes–Merton derivation is that it shows how to:**

* **A.** Increase the volatility of a portfolio.
* **B.** Find the average price of an asset.
* **C.** **Cancel out the random $dW_t$ term in a hedged portfolio to arrive at a deterministic PDE**. (**Correct**)
* **D.** Directly solve the Euler–Maruyama equation.

**2. Which mathematical term in the Black–Scholes PDE is directly related to the Itō Correction?**

* **A.** The $r\,f$ term.
* **B.** The $\frac{\partial f}{\partial t}$ term.
* **C.** The $\frac{\partial f}{\partial S}$ term.
* **D.** **The $\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 f}{\partial S^2}$ term (the second derivative with respect to price)**. (**Correct**)

---

#### Interview-Style Question

**Question:** In one sentence, summarize the central importance of Itō's Lemma to the field of quantitative finance.

**Answer Strategy:** Itō's Lemma is the essential stochastic chain rule that allows us to correctly model how the **value of a financial derivative** (which is a function of the stock price) changes over time, acknowledging that **volatility** contributes a **predictable, deterministic drift** that must be included in the derivative's valuation.

---

## 💡 Hands-On Simulation Projects (Chapter Conclusion) 🛠️

These projects are designed to implement the core stochastic calculus concepts, from the Wiener Process to the EM solution.

### Project 1: Simulating and Testing the Wiener Process

* **Goal:** Numerically verify the key properties of the Wiener Process.
* **Setup:** Choose $T=1.0$ and $N=10,000$ steps ($\Delta t = 10^{-4}$).
* **Steps:**
    1.  Generate the random increments $dW = \sqrt{\Delta t} Z$ where $Z \sim \mathcal{N}(0, 1)$.
    2.  Calculate the Wiener path $W_t$ by taking the cumulative sum of $dW$.
    3.  Verify two properties numerically:
        * The mean of $dW$ is approximately 0.
        * The variance of the final value $W_T$ is approximately $T=1.0$.
* ***Goal***: Confirm the core statistical properties of the driving noise source.

### Project 2: Visualizing the Order of Convergence (Strong)

* **Goal:** Visually demonstrate the strong convergence order $O(\sqrt{\Delta t})$ of the Euler–Maruyama method.
* **Setup:** Use the simple SDE $dS_t = \sigma dW_t$ ($\mu=0$) from $S_0=1.0$ to $T=1.0$ (known exact solution: $S_T = 1.0 + \sigma W_T$).
* **Steps:**
    1.  Choose a final noise value: $Z_{\text{final}}$ (used for the exact solution).
    2.  Run the EM simulation for $N=[10, 100, 1000, 10000]$ steps, ensuring the *total* accumulated noise is equal to $\sqrt{T} Z_{\text{final}}$ for all runs (a more complex strong-convergence requirement).
    3.  Calculate the **Absolute Error** $|S_T^{\text{exact}} - S_T^{\text{EM}}|$ for each $N$.
* ***Goal***: Plot the error versus $1/N$ ($\propto \Delta t$) on a log-log plot. The slope of the plot should be close to $0.5$ (since error $\propto \Delta t^{1/2}$), confirming the strong order of convergence.

### Project 3: The Itō Correction in Action (Numerical Check)

* **Goal:** Numerically confirm the presence of the $\mu - \frac{1}{2}\sigma^2$ drift in the log-price of GBM.
* **Setup:** Use GBM parameters: $S_0=100, \mu=0.10, \sigma=0.30, T=1.0$.
* **Steps:**
    1.  Run $M=10,000$ Euler–Maruyama simulations (or the exact GBM formula) and record the terminal price $S_T$ of each path.
    2.  Calculate the ensemble average of the final **log-price**: $\langle \ln(S_T / S_0) \rangle$.
    3.  Compare this numerical average with the theoretical expected log-return (the drift term in the Itō solution): $\mathbb{E}[\ln(S_T/S_0)] = (\mu - \tfrac{1}{2}\sigma^2)T$.
* ***Goal***: Show that the numerical average matches the term $(\mu - \tfrac{1}{2}\sigma^2)T$ (which is $0.055$) and *not* the simple expected return $\mu T$ (which is $0.10$), providing numerical evidence for the Itō correction.

### Project 4: Comparing EM and Exact GBM Solvers

* **Goal:** Compare the EM scheme (approximation) against the Exact formula for GBM, emphasizing the need for smaller $\Delta t$ in the EM solution.
* **Setup:** Use GBM parameters $\mu=0.10, \sigma=0.30, S_0=100, T=1.0$.
* **Steps:**
    1.  Run two simulations: **Simulation A** using the Exact GBM formula (one step, $N=1$) and **Simulation B** using the Euler–Maruyama formula (one step, $N=1$).
    2.  Run both A and B $M=10,000$ times and compare the mean terminal price $\langle S_T \rangle$ (which should both equal $S_0 e^{\mu T}$).
    3.  Now, run Simulation B again with $N=100$ steps ($\Delta t = 0.01$).
* ***Goal***: Demonstrate that the terminal mean of the one-step EM (Simulation B, $N=1$) will be slightly **biased** and not match the theoretical mean as accurately as the Exact formula, highlighting the discretization error inherent in EM for large $\Delta t$.
