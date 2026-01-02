# **Chapter 9: Black-Scholes-Merton Equation**

---

## **Introduction**

While stochastic differential equations (Chapter 8) correctly model the continuous randomness of asset price trajectories through geometric Brownian motion $dS_t = \mu S_t dt + \sigma S_t dW_t$, they leave a profound puzzle for financial valuation: how can we assign a **deterministic price** to derivatives whose payoffs depend on inherently random future prices? The naive approach—computing expected payoffs using Monte Carlo simulation—requires specifying the stock's expected return $\mu$, which varies across investors with different risk preferences. The 1973 breakthrough by Fischer Black, Myron Scholes, and Robert Merton resolved this paradox through a stunning insight: by constructing a **delta-hedged portfolio** that holds the derivative and dynamically shorts $\Delta = \frac{\partial V}{\partial S}$ units of the underlying stock, the random terms $dW_t$ driving both assets **cancel exactly**, creating a riskless portfolio whose return must equal the risk-free rate $r$ by no-arbitrage. This forces derivatives to satisfy a deterministic partial differential equation independent of $\mu$—the **Black–Scholes–Merton (BSM) equation**.

This chapter derives, analyzes, and numerically solves the BSM equation, revealing its deep connection to classical physics. We begin by applying Itō's Lemma to the option value $V(S,t)$ and constructing the delta-hedged portfolio $\Pi = V - \frac{\partial V}{\partial S} S$, showing how the cancellation of stochastic terms yields the BSM PDE: $\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0$. Through logarithmic transformation $x = \ln S$ and time reversal $\tau = T - t$, we expose the BSM equation's equivalence to the **heat diffusion equation** $\frac{\partial u}{\partial \tau} = \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2}$, where volatility acts as thermal diffusivity and option prices diffuse backward from the sharp payoff function at expiry. For European options, this yields the famous analytical solution $V = SN(d_1) - Ke^{-r(T-t)}N(d_2)$ via the heat kernel. American options, however, introduce a **free-boundary problem** requiring numerical solution.

By the end of this chapter, you will master the complete BSM framework: deriving the PDE through hedging arguments, transforming it to canonical diffusion form, and implementing the **Crank–Nicolson finite difference method** to solve American options with early exercise constraints $V(S,t) = \max(V_{\text{hold}}, V_{\text{intrinsic}})$. You will understand why the BSM equation eliminates dependence on expected returns (risk-neutral valuation), how the Itō correction term $\frac{1}{2}\sigma^2 S^2 V_{SS}$ emerges from quadratic variation, and why option pricing is fundamentally a diffusion problem. These techniques bridge stochastic modeling (Part I-II) to numerical PDE methods, preparing you for Chapter 10's exploration of reaction-diffusion systems in neuroscience, where similar mathematical structures govern action potential propagation.

---

## **Chapter Outline**

| **Sec.** | **Title** | **Core Ideas & Examples** |
|:---------|:----------|:--------------------------|
| **9.1** | From Random Paths to Deterministic Price | **Delta-hedged portfolio**: $\Pi = V - \Delta S$ with $\Delta = \frac{\partial V}{\partial S}$ cancels $dW_t$ terms. **No-arbitrage principle**: Riskless portfolio must earn $d\Pi = r\Pi dt$, eliminating dependence on $\mu$. **Risk-neutral valuation**: Option price independent of investor preferences. BSM PDE as parabolic diffusion equation. |
| **9.2** | The Derivation via Itō's Lemma | **Stochastic dynamics**: Stock $dS = \mu S dt + \sigma S dW_t$, option $dV = (\frac{\partial V}{\partial t} + \mu S V_S + \frac{1}{2}\sigma^2 S^2 V_{SS})dt + \sigma S V_S dW_t$. **Hedging construction**: $d\Pi = dV - \Delta dS$, choosing $\Delta = V_S$ eliminates random terms. **BSM PDE derivation**: $V_t + \frac{1}{2}\sigma^2 S^2 V_{SS} + rS V_S - rV = 0$, with $\mu$ term canceling. |
| **9.3** | BSM as Heat Equation | **Mathematical equivalence**: BSM PDE $\leftrightarrow$ heat diffusion equation. **Variable transformations**: Log-price $x = \ln S$, time-to-expiry $\tau = T-t$, exponential substitution $V = e^{\alpha x + \beta \tau}u(x,\tau)$ yields pure diffusion $u_\tau = \frac{1}{2}\sigma^2 u_{xx}$. **Analytical solution**: Heat kernel yields closed-form $V = SN(d_1) - Ke^{-r\tau}N(d_2)$ for European calls. |
| **9.4** | American Options with FDM | **Free-boundary problem**: Unknown optimal exercise boundary $S^*(t)$ analogous to Stefan problem. **Crank–Nicolson scheme**: Unconditionally stable, second-order accurate $\mathcal{O}(\Delta t^2, \Delta S^2)$, tridiagonal system $A\mathbf{V}^{n+1} = B\mathbf{V}^n$. **Early exercise constraint**: $V_i^n = \max(V_{\text{hold},i}^n, V_{\text{intrinsic}}(S_i))$ at every time step. Exercise frontier visualization. |
| **9.5** | Chapter Summary & Bridge | **Conceptual synthesis**: Stochastic GBM $\to$ Itō's Lemma $\to$ delta hedging $\to$ deterministic BSM PDE $\to$ diffusion analogy. **Risk neutralization**: Randomness eliminates itself through hedging, leaving only volatility and risk-free rate. Bridge to Chapter 10: From financial diffusion to neural reaction-diffusion (Hodgkin–Huxley model), nonlinear coupled ODEs for action potentials, transition from PDEs back to stiff ODE systems. |

---

## **9.1 From Random Paths to Deterministic Price**

-----

### **From Random Motion to Predictable Value**

In Chapter 8, we established the dynamics of a stock price $S_t$ using the **Geometric Brownian Motion (GBM)**, recognizing that prices follow a **random walk in continuous time**. The GBM SDE is expressed as:

$$dS_t = \mu S_t dt + \sigma S_t dW_t$$

Here, $\mu$ is the expected drift (return), $\sigma$ is the volatility, and $dW_t$ is the stochastic shock from the Wiener Process. This continuous randomness presents a challenge to **option pricing**, which requires estimating a value dependent on the uncertain future price.

While this expected value can be computed using computationally intensive statistical methods, such as **Monte Carlo simulation**, the 1973 breakthrough by **Fischer Black, Myron Scholes, and Robert Merton** provided a more profound insight: **option values obey a deterministic Partial Differential Equation (PDE)**.

-----

### **The Black–Scholes–Merton Breakthrough**

The core realization of the BSM model is that, despite the stochastic nature of the asset price, it is possible to **eliminate risk entirely** by constructing a dynamically managed portfolio. The model's key insight is that the randomness ($dW_t$ term) driving the stock price and the randomness driving the option value **cancel out** perfectly in a specifically constructed combination of the two assets.

This cancellation turns a probabilistic problem into a problem of determining a fair, **riskless** price in an arbitrage-free market.

-----

### **The Core Idea: The Delta-Hedged Portfolio**

To achieve this risk neutralization, a **delta-hedged portfolio ($\Pi$)** is constructed by holding one unit of the option (value $V$) and simultaneously shorting a specific number ($\Delta$) of the underlying stock ($S$):

$$\Pi = V - \Delta S$$

The goal is to choose the hedge ratio $\Delta$ such that the change in the portfolio value ($d\Pi$) over time $dt$ is purely deterministic (contains no $dW_t$ term).

The choice that perfectly cancels the stochastic terms is the option's **delta**:
$$\Delta = \frac{\partial V}{\partial S}$$

By setting $\Delta$ to this value, the resulting portfolio change, $d\Pi$, becomes deterministic. According to the **no-arbitrage principle**, any riskless portfolio must earn the **risk-free rate ($r$)**. This condition forces the deterministic change $d\Pi$ to equal the return on a risk-free bond, leading directly to the BSM PDE.

-----

### **The Black–Scholes–Merton PDE**

After rigorously applying **Itō’s Lemma** (Chapter 8) to the option value $V(S,t)$, substituting the result into $d\Pi = dV - \Delta dS$, and setting the $dW_t$ term to zero (by choosing $\Delta = \partial V/\partial S$), the following equation emerges:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + r S \frac{\partial V}{\partial S} - r V = 0$$

This is the **Black–Scholes–Merton (BSM) equation**. Notably, the PDE **does not depend on the stock's expected return ($\mu$)**. This independence is the core of **risk-neutral valuation**, stating that the option price is determined solely by volatility ($\sigma$), the risk-free rate ($r$), and the option's characteristics.

-----

### **The Physics Connection: Diffusion in Price Space**

The BSM equation is a **parabolic PDE**, which is mathematically identical to the well-known **heat or diffusion equation** from physics. This equivalence allows for the application of established computational methods:

| Concept | Physics (Heat Equation) | Finance (BSM Equation) |
| :--- | :--- | :--- |
| Field variable | Temperature $T(x,t)$ | Option price $V(S,t)$ |
| Spatial variable | Position $x$ | Asset price $S$ |
| Diffusion coefficient | Thermal diffusivity $\alpha$ | Volatility term $\frac{1}{2}\sigma^2 S^2$ |
| Initial condition | Temperature profile at $t=0$ | Payoff function at expiry $V(S,T)$ |

The option price, $V(S,t)$, behaves like a heat profile that **diffuses backward in financial time** from the sharp initial condition of the payoff function at expiration. This mathematical analogy allows for the reuse of numerical techniques like the **Finite Difference Method (FDM)** to solve the option pricing problem efficiently.

---

## **9.2 The Derivation: Applying Itō's Lemma**

We now formally derive the BSM PDE by applying **Itō's Lemma** to the option price function.

The Black–Scholes–Merton (BSM) equation is derived by combining the stochastic dynamics of the asset price (Chapter 8) with the powerful financial principle of **no arbitrage**. The objective is to construct a portfolio that is instantaneously **risk-free**, forcing its return to equal the risk-free rate ($r$).

-----

### **Step 1: The Stochastic Dynamics of the Asset and Option**

We begin with the stochastic differential equations governing the assets:

1.  **Stock Price Dynamics (Geometric Brownian Motion):** The price $S$ follows:
    $$dS = \mu S dt + \sigma S dW_t$$
    where $\mu$ is the expected return, $\sigma$ is the volatility, and $dW_t$ is the Wiener Process increment.

2.  **Option Value Dynamics (Itō's Lemma):** The derivative's value $V(S, t)$ is a function of the stochastic variable $S$. Applying **Itō's Lemma** (Chapter 8) to $V(S, t)$ yields its evolution:
    $$dV = \left(\frac{\partial V}{\partial t} + \mu S \frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}\right) dt + \sigma S \frac{\partial V}{\partial S} dW_t$$

!!! tip "Why Delta Hedging Eliminates Randomness"
    The key insight: when you hold $\Delta = \frac{\partial V}{\partial S}$ shares of stock and short one option, the instantaneous changes $dS_t$ and $dV_t$ **cancel out perfectly** in the combined portfolio. The $dW_t$ terms disappear, leaving only a deterministic drift that must equal the risk-free rate to prevent arbitrage.

-----

### **Step 2: Constructing the Delta-Hedged Portfolio**

To eliminate the risk associated with the random term ($dW_t$), we form a dynamically self-adjusting portfolio $\Pi$ that is **long one option** ($+V$) and **short a specific amount ($\Delta$) of the underlying stock** ($-S$):
$$\Pi = V - \Delta S$$

The change in the portfolio value ($d\Pi$) is:
$$d\Pi = dV - \Delta dS$$

Substituting the expressions for $dV$ and $dS$ into $d\Pi$ results in a long expression grouped by $dt$ (deterministic terms) and $dW_t$ (stochastic terms):

$$d\Pi = \underbrace{\left[\frac{\partial V}{\partial t} + \mu S \frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} - \Delta \mu S \right] dt}_{\text{Deterministic Terms}} + \underbrace{\left[\sigma S \frac{\partial V}{\partial S} - \Delta \sigma S \right] dW_t}_{\text{Stochastic Terms}}$$

-----

### **Step 3: Eliminating Randomness (The Hedge Ratio $\Delta$)**

The essence of the BSM breakthrough is choosing the hedge ratio $\Delta$ such that the stochastic term (the multiplier of $dW_t$) is exactly zero:

$$\sigma S \frac{\partial V}{\partial S} - \Delta \sigma S = 0$$

Solving for $\Delta$:
$$\Delta = \frac{\partial V}{\partial S}$$

This quantity $\Delta$ is the **option delta**, representing the sensitivity of the option price to a change in the stock price. By holding this precise ratio, the random market movement is neutralized, making $d\Pi$ purely deterministic.

Substituting $\Delta = \frac{\partial V}{\partial S}$ back into the $d\Pi$ equation results in the deterministic change in the risk-free portfolio:

$$d\Pi = \left[\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}\right] dt$$
*Note: The $\mu S \frac{\partial V}{\partial S}$ term, which contained the stock's expected return $\mu$, canceled with the $-\Delta \mu S$ term, which contained $\mu$ via the hedge ratio $\Delta$.*

-----

### **Step 4: The No-Arbitrage Condition**

Since the portfolio $\Pi$ is now riskless, the **no-arbitrage principle** dictates that its return must be the **risk-free rate ($r$)** over time $dt$:

$$d\Pi = r\Pi dt$$

Substituting the definition of $\Pi = V - \Delta S = V - S \frac{\partial V}{\partial S}$ into this condition:
$$d\Pi = r\left(V - S \frac{\partial V}{\partial S}\right) dt$$

-----

### **Step 5: Deriving the BSM PDE**

By equating the two deterministic expressions for $d\Pi$ (the Itō derivation and the no-arbitrage condition), we eliminate $dt$ and arrive at the final PDE:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} = rV - rS \frac{\partial V}{\partial S}$$

Rearranging to the standard form yields the **Black–Scholes–Merton Equation**:

$$\boxed{\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + r S \frac{\partial V}{\partial S} - r V = 0}$$

-----

### **The Profound Result: Risk-Neutral Valuation**

The final BSM PDE is notable because it is entirely independent of the stock's expected rate of return ($\mu$). This is the mathematical cornerstone of **risk-neutral valuation**: the price of the option does not depend on the market's collective forecast of future returns, only on the observable risk parameters ($\sigma$ and $r$).

---

## **9.3 The BSM Equation as a Heat Equation**

The **Black–Scholes–Merton (BSM) equation** is a parabolic partial differential equation that, despite its financial origin, is mathematically analogous to the classical **heat (diffusion) equation** from physics. This equivalence is crucial because it allows the problem of option pricing to be solved using established analytical and numerical techniques developed for diffusion processes.

-----

### **The Financial PDE and the Physics Analogy**

The BSM equation derived from the no-arbitrage principle is:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + r S \frac{\partial V}{\partial S} - r V = 0$$

The mathematical structure reveals its diffusion nature:

| Concept | Physics (Heat Equation) | Finance (BSM Equation) |
| :--- | :--- | :--- |
| **Field variable** | Temperature $T(x,t)$ | Option price $V(S,t)$ |
| **Spatial variable** | Position $x$ | Asset price $S$ |
| **Diffusion** | Thermal conductivity $\alpha$ | Volatility term $\frac{1}{2}\sigma^2 S^2$ |
| **Source/Sink term** | Heat loss (or gain) | Discounting term $-rV$ |

!!! example "The Heat Equation Analogy: Volatility as Thermal Diffusivity"
    Just as heat diffuses through a rod with diffusivity $\kappa$, option value diffuses through stock price space with diffusivity proportional to $\sigma^2$. High volatility = rapid diffusion of option value. The boundary conditions (strike price, expiration) act like temperature constraints at the rod's endpoints.

The analogy implies that uncertainty, driven by volatility ($\sigma$), causes the option value to **diffuse** or spread out over the space of possible stock prices ($S$) as time progresses.

-----

### **The Transformation to a Pure Diffusion Equation**

To fully expose this mathematical equivalence, the BSM equation is converted into the canonical form of the heat equation ($\frac{\partial u}{\partial \tau} = D \frac{\partial^2 u}{\partial x^2}$) using a series of variable substitutions.

-----

#### **Step 1: Logarithmic Price and Time Reversal**

Two initial transformations simplify the PDE by handling the multiplicative risk and reversing the time flow:
1.  **Logarithmic Price ($x$):** The substitution $x = \ln S$ linearizes the multiplicative terms ($\propto S \frac{\partial V}{\partial S}$ and $\propto S^2 \frac{\partial^2 V}{\partial S^2}$) by moving from price space to log-price space. This removes the variable coefficient $S^2$ from the second derivative term.
2.  **Time to Expiry ($\tau$):** The substitution $\tau = T - t$ flips the time axis such that the equation is solved forward in time-to-expiry ($\tau$), starting from the known payoff at $\tau=0$ (expiry, $t=T$).

Applying these changes transforms the BSM equation into a simpler, but not yet pure, diffusion form:

$$\frac{\partial V}{\partial \tau} = \frac{1}{2}\sigma^2 \frac{\partial^2 V}{\partial x^2} + \left(r - \frac{1}{2}\sigma^2\right)\frac{\partial V}{\partial x} - rV$$

-----

#### **Step 2: Removing the First Derivative and Decay Terms**

The remaining terms—the first derivative term $\propto \frac{\partial V}{\partial x}$ and the decay term $\propto -rV$—prevent the equation from being a standard heat equation.

A second, more complex substitution is used to eliminate these terms:
$$V(x, \tau) = e^{\alpha x + \beta \tau} u(x, \tau)$$

By carefully selecting constants $\alpha$ and $\beta$ to ensure the coefficients of $u_x$ and $u$ vanish, the equation reduces to the final, pure **heat (diffusion) equation** in the scaled variable $u$:

$$\frac{\partial u}{\partial \tau} = \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2}$$

-----

### **Initial and Boundary Conditions**

The option pricing problem is defined by its boundary and final conditions. When solving the PDE backward in time ($t \to 0$ from $t=T$), the **initial condition** in the transformed time $\tau$ is the known final payoff at expiry:

* **Initial Condition ($\tau=0$):** This is the option's payoff function at maturity. For a call option with strike $K$, this is $V(S, T) = \max(S - K, 0)$. In the heat analogy, this acts as a localized "heat pulse" that dissipates backward in time.
* **Boundary Conditions:** These ensure the solution behaves correctly at extreme prices ($S=0$ and $S \to \infty$). For example, a call option is worthless when the stock price is zero ($V(0,t)=0$).

-----

### **The Analytical Solution and Numerical Implications**

The pure diffusion form $u_\tau = \frac{1}{2}\sigma^2 u_{xx}$ has a known analytical solution involving the **heat kernel**. Transforming this solution back to the original variables ($V(S,t)$) yields the famous **Black–Scholes closed-form solution**:

$$V(S,t) = S N(d_1) - K e^{-r(T-t)} N(d_2)$$

This analytical solution relies on the **Gaussian weighting** function, which is mathematically equivalent to the diffusing heat kernel.

For cases where analytical solutions are impossible (such as for **American options**, due to the early exercise constraint), the established equivalence means that robust numerical methods like the **Finite Difference Method (FDM)**, including the Crank–Nicolson scheme, can be directly applied to the BSM PDE.

---

## **9.4 Simulation: Solving American Options with the Finite Difference Method (FDM)**

-----

### **The Free-Boundary Problem of American Options**

While European-style options, exercisable only at maturity, yield to the analytical closed-form solution of the Black–Scholes–Merton (BSM) equation, **American-style options** cannot be solved analytically. The difficulty arises because the holder has the right to exercise the option **at any time** before expiry ($t \le T$).

??? question "American vs European Options: When Does Early Exercise Matter?"
    For American **call options** on non-dividend-paying stocks, early exercise is **never optimal** (can be proven mathematically). But for American **put options**, early exercise becomes optimal when the stock price falls sufficiently below the strike—the option holder can capture the intrinsic value immediately rather than waiting and risking the stock price recovery. This asymmetry makes American puts significantly more complex to price than their European counterparts.

This early exercise feature introduces a **free-boundary problem** into the governing PDE:

1.  **The Governing Region:** The BSM PDE holds only in the region where it is optimal to **hold** the option.
2.  **The Exercise Region:** In the region where the intrinsic value (immediate payoff) is greater than the holding value, it is optimal to **exercise** immediately.
3.  **The Free Boundary ($\boldsymbol{S^*(t)}$):** The boundary separating these two regions—the optimal exercise price $S^* (t)$—is unknown and must be determined as part of the solution.

This situation is mathematically analogous to the **Stefan problem** in heat transfer, which models a moving phase front (like a melting or freezing boundary) in a diffusion system.

-----

### **The Inequality Constraint**

The value of an American option $V(S, t)$ must satisfy the standard BSM PDE in the holding region, but everywhere it must also satisfy the fundamental constraint that its value is never less than its immediate exercise payoff ($V_{\text{intrinsic}}$):

$$V(S, t) \ge V_{\text{intrinsic}}(S, t)$$

This makes the BSM equation an inequality known as a **Linear Complementarity Problem (LCP)**. The intrinsic value is the immediate payoff: $\max(S - K, 0)$ for a call or $\max(K - S, 0)$ for a put.

-----

### **Finite Difference Method (FDM) Framework**

Due to the free-boundary problem, numerical methods, particularly the **Finite Difference Method (FDM)**, are required. FDM discretizes the continuous space of the price ($S$) and time ($t$) into a finite grid.

The **Crank–Nicolson scheme** is the preferred method for solving the BSM PDE numerically. This semi-implicit scheme balances **stability** (unconditionally stable) and **accuracy** ($\mathcal{O}(\Delta t^2, \Delta S^2)$), making it ideal for the diffusion-like nature of the problem.

The solution process involves solving the BSM PDE **backward in time** from the known payoff at expiry ($t=T$ or $\tau=0$):

1.  **Discretization:** The PDE is approximated by a tridiagonal system of linear equations ($A \mathbf{V}^{n+1} = B \mathbf{V}^{n}$) at each time step.
2.  **Boundary Conditions:** Appropriate conditions (e.g., $V(0,t)=0$ for a call) are enforced at the edges of the price grid.

In an **implicit scheme** like **Crank–Nicolson**, the values of $V$ at time level $n+1$ are coupled, forming a system of linear equations:

$$
\mathbf{A} \mathbf{V}^{n+1} = \mathbf{B} \mathbf{V}^n
$$

Here is the algorithm structure:

```python
def crank_nicolson_american_option(S_grid, t_grid, r, sigma, K, option_type='put'):
    """
    Solve American option pricing using Crank-Nicolson FDM with early exercise.
    
    Parameters:
    - S_grid: Stock price grid (spatial discretization)
    - t_grid: Time grid (temporal discretization, backward from T to 0)
    - r: Risk-free rate
    - sigma: Volatility
    - K: Strike price
    - option_type: 'put' or 'call'
    """
    M = len(S_grid) - 1  # Number of spatial steps
    N = len(t_grid) - 1  # Number of time steps
    dS = S_grid[1] - S_grid[0]
    dt = t_grid[1] - t_grid[0]
    
    # Initialize option value at maturity (terminal condition)
    if option_type == 'put':
        V = np.maximum(K - S_grid, 0)  # Intrinsic value at expiration
    else:
        V = np.maximum(S_grid - K, 0)
    
    # Build tridiagonal matrices A and B for Crank-Nicolson
    alpha = 0.25 * dt * (sigma**2 * (np.arange(M+1)**2) - r * np.arange(M+1))
    beta = -0.5 * dt * (sigma**2 * (np.arange(M+1)**2) + r)
    gamma = 0.25 * dt * (sigma**2 * (np.arange(M+1)**2) + r * np.arange(M+1))
    
    # March backward in time
    for n in range(N):
        # Solve A * V_new = B * V_old (tridiagonal system)
        V_new = solve_tridiagonal_system(A, B, V)
        
        # Apply early exercise constraint (American option)
        if option_type == 'put':
            V_intrinsic = np.maximum(K - S_grid, 0)
        else:
            V_intrinsic = np.maximum(S_grid - K, 0)
        
        V = np.maximum(V_new, V_intrinsic)  # Choose max of hold vs exercise
    
    return V  # Option value at t=0 for all stock prices
```

-----

### **Implementing the Early Exercise Constraint**

The crucial step in the FDM algorithm for American options is the enforcement of the early exercise constraint **at every single time step**. After the Crank–Nicolson scheme is solved (which yields the value assuming the option is held), the solution must be immediately checked against the intrinsic value:

$$V_i^{n} = \max\left(V_{\text{hold}, i}^{n}, V_{\text{intrinsic}}(S_i, t_n)\right)$$

where $V_{\text{hold}}$ is the value computed by the PDE discretization. This selection ensures that the calculated price $V$ is always the maximum of the continuation value (holding) and the exercise value (intrinsic).

This constant application of the $\max$ function dynamically captures the optimal exercise strategy and correctly tracks the location of the **free boundary** $S^*(t)$ as the solution propagates backward from expiry to the present.

-----

### **Financial Interpretation of the Solution**

The FDM simulation provides a complete picture of the option's value surface $V(S,t)$ over the entire domain of price and time.

* **Early Exercise Premium:** The numerical American price will always be greater than or equal to the corresponding European price ($V_{\text{American}} \ge V_{\text{European}}$) due to the added flexibility, with the difference representing the **early exercise premium**.
* **Exercise Frontier:** Visualizing the boundary where $V(S,t) = V_{\text{intrinsic}}(S,t)$ reveals the **optimal exercise strategy**—the price at which a rational investor should exercise the option immediately. For American put options, this boundary is active when the stock price is low, signaling that the option should be exercised before discounting erodes the large intrinsic value.

This numerical solution demonstrates how the combination of a continuous stochastic model (GBM), a deterministic parabolic PDE (BSM), and a discrete boundary constraint (FDM) provides a complete model for complex derivative valuation.

---

## **9.5 Chapter Summary & Bridge to Chapter 10**

Chapter 9 completed the derivation and computational solution of the **Black–Scholes–Merton (BSM) equation**, marking a synthesis of stochastic mathematics and deterministic partial differential equations (PDEs). The central achievement was transforming the inherent randomness of asset prices into a predictable valuation framework.

```mermaid
flowchart TD
    A[Stock Price: GBM Model<br/>dS = μS dt + σS dW] --> B[Apply Itō's Lemma<br/>to Option Value V(S,t)]
    B --> C[Construct Delta-Hedged Portfolio<br/>Π = V - ΔS]
    C --> D[Compute Portfolio Change dΠ]
    D --> E{Cancel dW Terms<br/>via Δ = ∂V/∂S}
    E --> F[Deterministic Portfolio Drift<br/>No Randomness Remaining]
    F --> G[Apply No-Arbitrage Condition<br/>dΠ = rΠ dt]
    G --> H[Black-Scholes-Merton PDE<br/>∂V/∂t + ½σ²S²∂²V/∂S² + rS∂V/∂S - rV = 0]
    H --> I{Option Type?}
    I -->|European| J[Analytical Solution<br/>BSM Formula]
    I -->|American| K[Numerical Solution<br/>Crank-Nicolson FDM]
    K --> L[Free-Boundary Problem<br/>Early Exercise Constraint]
    
    style H fill:#e1f5ff
    style J fill:#d4edda
    style K fill:#fff3cd
```

-----

### **The Unification of Stochastic and Deterministic Dynamics**

The BSM framework hinges on the dynamic interplay of several mathematical components:

* **Stochastic Foundation:** The asset price dynamics, modeled as **Geometric Brownian Motion** ($dS = \mu S dt + \sigma S dW_t$), required **Itō's Lemma** to correctly calculate the change in the option's value ($dV$).
* **Risk Neutralization:** By constructing a **delta-hedged portfolio** ($\Pi = V - \frac{\partial V}{\partial S} S$), the random $dW_t$ term in the dynamics of $\Pi$ was perfectly **canceled out**.
* **Deterministic Law:** The no-arbitrage condition forced the resulting riskless change ($d\Pi$) to equal the risk-free return ($r\Pi dt$), leading to the final BSM PDE:
    $$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + r S \frac{\partial V}{\partial S} - r V = 0.$$
* **Physical Equivalence:** This PDE is mathematically identical to the **Heat (Diffusion) Equation**, confirming that option pricing is fundamentally a problem of **diffusion in price space**, where **volatility ($\sigma$) acts as the thermal diffusivity**.

-----

### **Solving the Free Boundary Problem**

The BSM equivalence to the diffusion equation allowed the application of robust numerical solvers:

* **European Options** were solvable via the analytical closed form, derived from the fundamental solution (the heat kernel).
* **American Options** presented a complex **free-boundary problem** due to the optimal early exercise constraint. This was handled numerically by solving the BSM PDE backward in time using the stable **Crank–Nicolson FDM scheme**, while enforcing the inequality constraint: $V(S, t) = \max\big(V_{\text{hold}}(S, t), V_{\text{intrinsic}}(S, t)\big)$ at every time step.

This numerical technique for the American option is analogous to tracking a **moving phase boundary** (e.g., a melting front) in physical diffusion problems.

-----

### **Bridge to Chapter 10: From Financial Diffusion to Neural Dynamics**

The preceding chapters demonstrated that systems governed by **diffusion** and **non-linear feedback** can be found across different disciplines.

| Domain | Dynamic Variable | Governing Process | Key Mathematical Form |
| :--- | :--- | :--- | :--- |
| **Finance** | Option Price ($V$) | Diffusion under Volatility | BSM PDE |
| **Neuroscience** | Membrane Voltage ($V_m$) | Reaction–Diffusion of Ions | Coupled Nonlinear ODEs |

In the financial model, **information** diffuses through the market as price fluctuations. In the upcoming **Biological model**, **ions** diffuse across a neural membrane, generating electrical impulses.

**Chapter 10** will introduce the **Hodgkin–Huxley model**, which describes how a neuron generates an **action potential** (a nerve impulse). This system shifts the computational focus back from **PDEs** to a highly nonlinear system of **coupled Ordinary Differential Equations (ODEs)**:

$$C_m \frac{dV_m}{dt} = -I_{\text{ion}} + I_{\text{stim}}$$

The simulation challenge here is no longer diffusion across space but the intricate, non-linear time-evolution of **voltage ($V_m$)** and **gating variables ($m, h, n$)** that control ion channels, requiring specialized numerical integration techniques for stiff ODEs. The continuous narrative remains: local interactions (ion current) lead to emergent global phenomena (neural signal).

These techniques will continue to serve as mathematical foundations for understanding **complex, high-dimensional stochastic systems** encountered in the subsequent chapters.

---

## **References**

1. **Black, F., & Scholes, M. (1973).** "The Pricing of Options and Corporate Liabilities." *Journal of Political Economy*, 81(3), 637–654. — The original groundbreaking paper deriving the BSM formula.

2. **Merton, R. C. (1973).** "Theory of Rational Option Pricing." *Bell Journal of Economics and Management Science*, 4(1), 141–183. — Extends Black-Scholes framework with rigorous no-arbitrage theory.

3. **Hull, J. C. (2022).** *Options, Futures, and Other Derivatives* (11th ed.). Pearson. — Comprehensive textbook covering option pricing theory and numerical methods.

4. **Wilmott, P. (2006).** *Paul Wilmott on Quantitative Finance* (2nd ed.). Wiley. — Detailed treatment of BSM derivation, Greeks, and numerical PDE methods for option pricing.

5. **Crank, J., & Nicolson, P. (1947).** "A Practical Method for Numerical Evaluation of Solutions of Partial Differential Equations of the Heat-Conduction Type." *Proceedings of the Cambridge Philosophical Society*, 43(1), 50–67. — Original paper on the Crank-Nicolson implicit scheme.

6. **Brennan, M. J., & Schwartz, E. S. (1977).** "The Valuation of American Put Options." *Journal of Finance*, 32(2), 449–462. — Early work on numerical methods for American option pricing with free-boundary problems.

7. **Shreve, S. E. (2004).** *Stochastic Calculus for Finance II: Continuous-Time Models*. Springer. — Rigorous treatment of martingale pricing theory and BSM derivation using risk-neutral valuation.

8. **Tavella, D., & Randall, C. (2000).** *Pricing Financial Instruments: The Finite Difference Method*. Wiley. — Practical guide to implementing FDM for option pricing, including American options.

9. **Duffy, D. J. (2013).** *Finite Difference Methods in Financial Engineering: A Partial Differential Equation Approach*. Wiley. — Comprehensive coverage of numerical PDE techniques for computational finance.

10. **Glasserman, P. (2003).** *Monte Carlo Methods in Financial Engineering*. Springer. — Comparison of PDE vs Monte Carlo approaches for option pricing; highlights when each method excels.

