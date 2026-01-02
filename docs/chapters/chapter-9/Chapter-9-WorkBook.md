## 💰 Chapter 9: Finance III: Black-Scholes-Merton (BSM) (Workbook)

The goal of this chapter is to connect stochastic calculus (Chapter 8) with deterministic PDE solvers (Volume I) to derive and numerically solve the **Black–Scholes–Merton (BSM) equation**, which prices derivatives by eliminating market risk.

| Section | Topic Summary |
| :--- | :--- |
| **9.1** | Chapter Opener: From Random Paths to Deterministic Price |
| **9.2** | The Derivation: Applying Itō’s Lemma |
| **9.3** | The BSM Equation as a Heat Equation |
| **9.4** | Simulation: Solving American Options with FDM |
| **9.5** | Chapter Summary & Bridge to Chapter 10 |

***

### 9.1 From Random Paths to Deterministic Price

> **Summary:** The **Black–Scholes–Merton (BSM) breakthrough** showed that despite the stochastic nature of asset prices ($dS_t = \mu S_t\,dt + \sigma S_t\,dW_t$), option values obey a **deterministic Partial Differential Equation (PDE)**. This is achieved by constructing a **delta-hedged portfolio** where the random components ($dW_t$) exactly **cancel**.

#### Section Detail

The delta-hedged portfolio ($\Pi = V - \Delta S$) is designed to be **riskless**. By the **no-arbitrage principle**, this riskless portfolio must earn the **risk-free rate** ($r$). This conversion of a stochastic problem into a deterministic condition is the core of the BSM model. The resulting PDE is **parabolic**, mathematically equivalent to the Heat/Diffusion Equation.

#### Quiz Questions

**1. The breakthrough insight of the BSM model that allows for deterministic pricing is the discovery that:**

* **A.** Volatility is always zero.
* **B.** The expected return ($\mu$) is always equal to the risk-free rate ($r$).
* **C.** **The random component ($dW_t$) of the stock and derivative dynamics can be made to cancel in a hedged portfolio**. (**Correct**)
* **D.** All stock prices follow a normal distribution.

**2. The resulting BSM equation is classified as a parabolic PDE, which is mathematically equivalent to which fundamental equation from physics?**

* **A.** The Wave Equation.
* **B.** The Navier-Stokes Equation.
* **C.** The **Heat/Diffusion Equation**. (**Correct**)
* **D.** The Schrödinger Equation.

---

#### Interview-Style Question

**Question:** Briefly define the concept of a **delta-hedged portfolio** in the BSM context and explain its primary purpose.

**Answer Strategy:** A delta-hedged portfolio is a combination of holding one unit of the derivative ($V$) and simultaneously shorting $\Delta = \partial V / \partial S$ units of the underlying asset ($S$). Its primary purpose is to **neutralize (cancel out)** the random $dW_t$ term in the portfolio's change in value, thereby creating a riskless asset. The subsequent condition that this riskless asset must earn the risk-free rate ($d\Pi = r\Pi\,dt$) leads directly to the BSM PDE.

---

***

### 9.2 The Derivation: Applying Itō’s Lemma

> **Summary:** The BSM equation is derived by applying **Itō’s Lemma** to the derivative value $V(S, t)$, substituting the result into the change in the portfolio $d\Pi = dV - \Delta\,dS$, and setting the random part to zero by choosing $\Delta = \partial V / \partial S$. The resulting deterministic equation, combined with the **no-arbitrage condition** ($d\Pi = r\Pi\,dt$), produces the BSM PDE, which notably **does not depend on the stock's expected return ($\mu$)**.

#### Section Detail

The cancellation of the random term is perfect, leaving $d\Pi$ purely deterministic. The profound result is that the final price depends only on $\sigma$ (volatility), $r$ (risk-free rate), $S$ (price), $K$ (strike), and $T$ (time to expiry), but not $\mu$. This is why the BSM model works under the risk-neutral measure.

#### Quiz Questions

**1. The value chosen for the hedge ratio $\Delta$ (the number of shares to hold against one option) that ensures the portfolio is risk-free is:**

* **A.** $\Delta = \sigma / r$.
* **B.** $\Delta = r / \sigma$.
* **C.** $\Delta = \frac{\partial V}{\partial t}$.
* **D.** $\Delta = \frac{\partial V}{\partial S}$ (the option's delta). (**Correct**)

**2. A profound observation of the final BSM PDE is that the price $V(S,t)$ is independent of which variable?**

* **A.** Volatility ($\sigma$).
* **B.** The risk-free rate ($r$).
* **C.** **The stock's expected rate of return ($\mu$)**. (**Correct**)
* **D.** The strike price ($K$).

---

#### Interview-Style Question

**Question:** The deterministic term in the portfolio change $d\Pi$ still contains the $\mu S \partial V / \partial S$ term from Itō's Lemma. Explain how $\mu$ is ultimately eliminated from the final BSM PDE.

**Answer Strategy:** The $\mu$ term is eliminated in the final step where the **no-arbitrage condition** is imposed. The final step equates the deterministic change derived from Itō's Lemma with the return of a riskless bond:
$$d\Pi_{\text{Ito}} = d\Pi_{\text{Risk-Free}}$$
$$
\left(\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}\right) dt = r(V - \Delta S) dt
$$
By substituting $\Delta = \partial V / \partial S$, the stock's expected return $\mu$ is no longer necessary to describe the derivative's value, as its effects are perfectly offset by the hedging strategy, and the portfolio must simply return the risk-free rate $r$.

---

***

### 9.3 The BSM Equation as a Heat Equation

> **Summary:** The BSM equation is mathematically equivalent to the standard **Heat (Diffusion) Equation**. This equivalence is revealed through a **change of variables** ($x=\ln S$, $\tau=T-t$) and a further function substitution ($V = e^{\alpha x + \beta \tau} u$), which simplifies the BSM PDE to the pure diffusion form $u_{\tau} = \frac{1}{2}\sigma^2 u_{xx}$.

#### Section Detail

The BSM equation is a parabolic PDE where volatility ($\frac{1}{2}\sigma^2 S^2$) acts as the diffusion coefficient. The transformation allows the problem to be solved either analytically (resulting in the famous BSM closed-form solution) or numerically using standard FDM methods. The initial condition for the solution is the option's payoff at expiry ($V(S,T)$), which acts as a localized **"heat pulse"** that diffuses backward in financial time.

#### Quiz Questions

**1. Which change of variables is necessary to transform the BSM equation into a pure diffusion equation?**

* **A.** Substituting $\sigma$ with $\mu$ and $r$ with $\sigma$.
* **B.** **Transforming the price $S$ into a log-space coordinate $x = \ln S$ and time $t$ into time-to-expiry $\tau = T - t$**. (**Correct**)
* **C.** Multiplying the entire equation by the time step $\Delta t$.
* **D.** Setting the diffusion term to zero.

**2. In the Heat Equation analogy, the option value $V(S,t)$ corresponds to the physical quantity of:**

* **A.** Heat flux.
* **B.** Mass density.
* **C.** **Temperature $T(x,t)$**. (**Correct**)
* **D.** Thermal conductivity.

---

#### Interview-Style Question

**Question:** The transformation of the BSM equation often involves two major steps: $x=\ln S$ and $V = e^{\alpha x + \beta \tau} u$. Explain the purpose of the second, more complex substitution ($V = e^{\alpha x + \beta \tau} u$).

**Answer Strategy:** The first substitution ($x=\ln S$) linearizes the multiplicative randomness, but the transformed PDE still contains a **first derivative term ($V_x$)** and a **decay term ($-rV$)**. The substitution $V = e^{\alpha x + \beta \tau} u$ is specifically chosen to eliminate these remaining unwanted terms. By correctly choosing constants $\alpha$ and $\beta$, the resulting equation simplifies to the classic, pure diffusion form ($u_{\tau} = \frac{1}{2}\sigma^2 u_{xx}$), allowing the use of standard heat equation solution techniques.

---

***

### 9.4 Simulation: Solving American Options with FDM

> **Summary:** Analytical solutions for the BSM PDE exist only for European options. **American options** (exercisable anytime) require **numerical methods** like the **Finite Difference Method (FDM)** because they present a complex **free-boundary problem**. The numerical solution requires solving the Crank–Nicolson discretized PDE backward in time while explicitly applying the **early exercise constraint** at every time step.

#### Section Detail

The FDM discretizes the BSM PDE onto a grid of price ($S$) and time ($t$). The stable and accurate **Crank–Nicolson scheme** is used to solve the resulting tridiagonal matrix system backward in time. The critical step for American options is the **constraint enforcement**: at every point ($S_i, t_n$), the calculated option value $V$ must be the maximum of the holding value (from the PDE solution) and the intrinsic value ($\max(S-K, 0)$ or $\max(K-S, 0)$).

#### Quiz Questions

**1. The primary feature that makes the American option pricing problem intractable for analytical solutions (and thus necessitates FDM) is that it is a:**

* **A.** Log-normal distribution problem.
* **B.** $\mathcal{O}(N^2)$ complexity problem.
* **C.** **Free-boundary problem**. (**Correct**)
* **D.** Time-reversible problem.

**2. Which widely-used FDM scheme is typically favored for solving the BSM PDE due to its balance of stability and $\mathcal{O}(\Delta t^2, \Delta S^2)$ accuracy?**

* **A.** The Explicit Euler scheme.
* **B.** The Implicit Euler scheme.
* **C.** The **Crank–Nicolson scheme**. (**Correct**)
* **D.** The Milstein scheme.

---

#### Interview-Style Question

**Question:** In the FDM algorithm for American options, why is the optimal exercise boundary problem solved by taking the **maximum** of the calculated PDE value and the option’s intrinsic value?

**Answer Strategy:** The FDM solves the PDE for the option's value **if it is held** (the time value). However, the holder has the right to exercise, which yields the intrinsic value. At any given time, a rational holder will choose the path that yields the highest value. Therefore, the option's true value must be $\max(V_{\text{hold}}, V_{\text{intrinsic}})$. The point where this maximum switches defines the optimal exercise frontier (the free boundary).

---

***

### 9.5 Chapter Summary & Bridge to Chapter 10

> **Summary:** Chapter 9 demonstrated the power of mathematics to transform a **stochastic problem** into a **deterministic PDE** via the BSM framework. The BSM equation’s equivalence to the **diffusion equation** confirms that option pricing is a problem of **diffusion in price space**. Numerical solutions (FDM) are essential for handling complex constraints like the **free boundary of American options**.

#### Section Detail

The synthesis of Itō calculus, the no-arbitrage principle, and PDE methods provides a complete model for derivative valuation. The challenge of the American option (solving a PDE with a moving boundary) is analogous to complex problems in physics (like the Stefan problem), reinforcing the interdisciplinary nature of computational methods. The next chapter shifts from financial diffusion to ion diffusion in neurons.

#### Quiz Questions

**1. The analogy between the BSM equation and the Heat Equation means that market **volatility ($\sigma$)** in finance corresponds physically to which property?**

* **A.** Time (t).
* **B.** **Thermal conductivity (or diffusivity, $\alpha$)**. (**Correct**)
* **C.** Temperature (T).
* **D.** Energy (E).

**2. Which core concept, shared across finance, thermodynamics, and physics, is essential for transforming the stochastic stock price SDE into the deterministic BSM PDE?**

* **A.** The Law of Large Numbers.
* **B.** The Central Limit Theorem.
* **C.** The concept that **randomness can be neutralized in a self-consistent system**. (**Correct**)
* **D.** The need for a Monte Carlo simulation.

---

#### Interview-Style Question

**Question:** In the context of the BSM model, how is the mathematical concept of a **moving phase boundary** from physics relevant to financial decision-making?

**Answer Strategy:** The moving phase boundary is a direct analogy to the **optimal early exercise frontier** ($S^*(t)$) for an American option.
* **Phase Boundary (Physics):** Separates two states (e.g., solid/liquid) where the diffusion equation applies on one side and a constraint applies on the other.
* **Exercise Frontier (Finance):** Separates the region where it is optimal to **hold** the option (where the BSM PDE applies) from the region where it is optimal to **exercise** immediately (where the constraint $V = V_{\text{intrinsic}}$ applies). Numerically solving the American option requires tracking this moving boundary.

---

## 💡 Hands-On Simulation Projects (Chapter Conclusion) 🛠️

These projects are designed to implement the core BSM solution via numerical PDE methods.

### Project 1: Testing the Analytical BSM Solution

* **Goal:** Calculate the analytical BSM price for a European call and observe the effect of volatility.
* **Setup:** Use $S=100, K=100, r=0.05, T=1.0$.
* **Steps:**
    1.  Implement the BSM analytical formula (using $d_1, d_2$, and the standard normal CDF).
    2.  Calculate the price $V_{\text{call}}$ for $\sigma=0.10$ and $\sigma=0.50$.
    3.  Calculate the **option Vega** ($\frac{\partial V}{\partial \sigma}$), which is the sensitivity of price to volatility.
* ***Goal***: Show that the option price significantly increases with $\sigma$ and confirm that the Vega is always positive (an option is always worth more when the future is more uncertain).

### Project 2: Implementing the Forward Euler FDM Scheme for BSM

* **Goal:** Implement the simplest FDM scheme to solve the BSM PDE, solving backward in time.
* **Setup:** Use $S_{\text{max}} = 200, K=100, r=0.05, \sigma=0.20, T=1.0$. Discretize the grid $N_S=100, N_t=500$.
* **Steps:**
    1.  Implement the **Explicit (Forward Euler)** discretization scheme for the BSM PDE (which solves the grid $V^n$ using $V^{n+1}$ only).
    2.  Set the final payoff $V(S,T) = \max(S-K, 0)$ as the initial condition.
    3.  Iterate backward in time, $n = N_t \to 0$, applying the appropriate boundary conditions (e.g., $V(0,t)=0$).
* ***Goal***: Produce a numerical price $V(S_0, 0)$ and compare it against the analytical BSM price (Project 1) to check the numerical accuracy.

### Project 3: Stability Check for the Explicit FDM Scheme

* **Goal:** Demonstrate the numerical instability inherent in the Explicit FDM scheme.
* **Setup:** Use the same parameters as Project 2, but intentionally choose an unstable combination of grid parameters, e.g., $\Delta t$ too large relative to $\Delta S$. (Explicit FDM requires a stability condition: $\Delta t \le \Delta S^2 / (\sigma^2 S^2)$).
* **Steps:**
    1.  Run the Explicit FDM solver with the unstable parameters.
    2.  Monitor the calculated option values $V$.
* ***Goal***: Show that the $V$ values quickly "blow up" to non-physical, oscillating, or infinite values, illustrating the critical importance of scheme stability in PDE solvers.

### Project 4: Modeling the Early Exercise Constraint (American Put)

* **Goal:** Numerically enforce the early exercise constraint to price an American option.
* **Setup:** Use the stable Crank-Nicolson scheme (or a simple Implicit scheme) for discretization. $S_{\text{max}} = 200, K=100, r=0.05, \sigma=0.20, T=1.0$.
* **Steps:**
    1.  Implement the Crank-Nicolson discretization (or use a library function for the implicit solve).
    2.  After solving the system at each backward time step $n$, apply the early exercise check: $V_i^{n} = \max(V_i^{n}, K - S_i)$.
* ***Goal***: Compare the resulting American Put price $V_{\text{American}}$ with the corresponding European Put price $V_{\text{European}}$ (calculated without the early exercise constraint). Show that $V_{\text{American}} \ge V_{\text{European}}$ (the early exercise premium).


