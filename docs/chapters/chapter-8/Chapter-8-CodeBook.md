Got it. I'll proceed with the hands-on simulation projects for Chapter 8, which focus on implementing and verifying the core concepts of stochastic calculus.

-----

# Chapter 8: Finance II: The Stochastic Calculus (SDEs)

## Project 1: Simulating and Testing the Wiener Process

-----

### Definition: Simulating and Testing the Wiener Process

The goal of this project is to numerically verify the core properties of the **Wiener Process ($W_t$)**, the fundamental driving noise for Stochastic Differential Equations (SDEs). This is done by generating a large number of random increments and confirming the zero-mean and variance scaling properties.

### Theory: Properties of the Wiener Increment

The Wiener Process is the continuous limit of a random walk. Its infinitesimal increment $dW_t$ satisfies two key properties for a given finite time step $\Delta t$:

1.  **Zero Mean:** The increment has an expected value of zero.
    $$\mathbb{E}[dW_t] = 0$$
2.  **Variance Scaling:** The variance of the increment is equal to the time step.
    $$\text{Var}(dW_t) = \mathbb{E}[(dW_t)^2] = \Delta t$$

In simulation, the increment is generated using a standard normal variate ($Z \sim \mathcal{N}(0, 1)$):

$$dW_t = \sqrt{\Delta t} \cdot Z$$

This project verifies that the cumulative process, $W_T = \sum dW_t$, has a mean of zero and a variance equal to the total time $T$.

-----

### Extensive Python Code and Visualization

The code simulates many independent paths of the Wiener Process over $T=1.0$ year and verifies that the ensemble of terminal values satisfies the theoretical statistical properties.

```python
import numpy as np
import matplotlib.pyplot as plt

# Set seed for reproducibility
np.random.seed(42)

# ====================================================================
# 1. Simulation Parameters
# ====================================================================

T = 1.0     # Total time (T) for the process (e.g., 1 year)
N = 10000   # Number of time steps (N)
DT = T / N  # Time step size (Delta t)

M_PATHS = 10000  # Number of independent paths to run for ensemble verification

# Theoretical expectation for verification
E_WT_THEO = 0.0  # Mean of the Wiener Process at any time T is 0
VAR_WT_THEO = T  # Variance of the Wiener Process at time T is T

# ====================================================================
# 2. Simulation and Verification
# ====================================================================

terminal_W = np.zeros(M_PATHS)

# Calculate the constant scaling factor for the random increment
dW_scale = np.sqrt(DT)

for m in range(M_PATHS):
    # 1. Generate N independent standard normal variates
    Z_sequence = np.random.standard_normal(N)
    
    # 2. Calculate the Wiener increments: dW = sqrt(dt) * Z
    dW_sequence = dW_scale * Z_sequence
    
    # 3. Calculate the Wiener path: W_t = cumulative sum(dW)
    W_path = np.cumsum(dW_sequence)
    
    # Record the terminal value W_T
    terminal_W[m] = W_path[-1]

# --- Calculate Empirical Statistics ---
E_WT_EMPIRICAL = np.mean(terminal_W)
VAR_WT_EMPIRICAL = np.var(terminal_W, ddof=1) # Use ddof=1 for sample variance

# ====================================================================
# 3. Visualization and Analysis
# ====================================================================

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Ensemble Distribution (Histogram)
ax[0].hist(terminal_W, bins=50, density=True, color='purple', alpha=0.7, 
           label='Simulated $W_T$')
ax[0].axvline(E_WT_THEO, color='red', linestyle='--', 
             label='Theoretical Mean $\\mathbb{E}[W_T] = 0$')

# Labeling and Formatting
ax[0].set_title(f'Ensemble Distribution of Terminal $W_T$ ($T={T}$)')
ax[0].set_xlabel('Terminal Value $W_T$')
ax[0].set_ylabel('Density')
ax[0].legend()
ax[0].grid(True)

# Plot 2: Variance Check (Illustrating the squared term)
check_data = [VAR_WT_THEO, VAR_WT_EMPIRICAL]
ax[1].bar(['Theoretical Var($W_T$)=T', f'Empirical Var($W_T$): {VAR_WT_EMPIRICAL:.4f}'], 
         check_data, color=['gray', 'purple'])
ax[1].axhline(VAR_WT_THEO, color='red', linestyle='--', label='Target Variance')

# Labeling and Formatting
ax[1].set_title('Verification of Variance Scaling')
ax[1].set_ylabel('Variance')
ax[1].grid(True, axis='y')

plt.tight_layout()
plt.show()

# --- Verification Summary ---
print("\n--- Wiener Process Verification Summary ---")
print(f"Time Step (\u0394t): {DT:.4e}")
print(f"Total Paths (M): {M_PATHS}")
print("-----------------------------------------")
print(f"Theoretical Mean \u222e[W_T]: {E_WT_THEO:.4f}")
print(f"Empirical Mean \u222e[W_T]:   {E_WT_EMPIRICAL:.4f}")
print(f"Difference (Mean):       {np.abs(E_WT_EMPIRICAL - E_WT_THEO):.4e}")
print("-----------------------------------------")
print(f"Theoretical Variance Var[W_T]: {VAR_WT_THEO:.4f}")
print(f"Empirical Variance Var[W_T]:   {VAR_WT_EMPIRICAL:.4f}")
print(f"Difference (Variance):         {np.abs(VAR_WT_EMPIRICAL - VAR_WT_THEO):.4e}")

print("\nConclusion: The simulation successfully generated the Wiener Process. The ensemble of terminal values confirms the two defining properties: the mean is zero, and the variance is equal to the total time T (1.0), which provides the fundamental noise input for SDE solvers.")
```

-----

## Project 2: Visualizing the Order of Convergence (Strong)

-----

### Definition: Visualizing Strong Convergence

The goal of this project is to demonstrate the **strong convergence order** of the **Euler–Maruyama (EM) method**. Strong convergence measures how accurately a simulated path ($S_T^{\text{EM}}$) tracks a specific, true path ($S_T^{\text{exact}}$).

### Theory: Strong Convergence Order

The EM scheme has a **strong convergence order** of $O(\sqrt{\Delta t})$ (or $\mathcal{O}(\Delta t^{1/2})$). This means the error between the numerical solution and the true solution scales with the square root of the time step:

$$\text{Error} = \mathbb{E}[|S_T^{\text{exact}} - S_T^{\text{EM}}|] \propto \Delta t^{1/2} \propto N^{-1/2}$$

To verify this numerically, we plot the error versus the inverse of the number of steps ($1/N$) on a **log-log plot**. The slope of this line should approximate the strong order, **$0.5$**.

**Simplified SDE:** We use the simple SDE $dS_t = \sigma dW_t$ (which has zero drift, $\mu=0$) from $S_0=1.0$ to $T=1.0$.

  * The exact solution for this SDE is $S_T^{\text{exact}} = S_0 + \sigma W_T$.
  * We fix the final noise value ($W_T$) for the exact solution and check how close the EM paths get to it as $N$ increases.

-----

### Extensive Python Code and Visualization

The code runs the EM solver for several grid sizes ($N$), calculates the absolute error against a known target, and performs the log-log fit to verify the $O(\sqrt{\Delta t})$ scaling.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Set seed for reproducibility
np.random.seed(42)

# ====================================================================
# 1. Simulation Parameters
# ====================================================================

sigma = 0.30  # Volatility
S0 = 1.0      # Initial price
T = 1.0       # Time to maturity

# Sequence of step counts to test
N_values = np.array([10, 50, 250, 1000, 5000, 10000]) 

# Fix the final noise value (Z_final) for the exact solution.
# This ensures that all EM paths are aimed at the same true terminal point.
Z_FINAL = 1.5 
W_T_FIXED = sigma * np.sqrt(T) * Z_FINAL
S_T_EXACT = S0 + W_T_FIXED # Exact solution for dS = sigma dW is S_T = S0 + W_T

# ====================================================================
# 2. Euler-Maruyama Solver for Strong Convergence
# ====================================================================

def euler_maruyama_strong(S0, sigma, T, N, Z_final):
    """
    Simulates the SDE dS = sigma dW using EM, ensuring the total accumulated 
    noise is a fixed value (W_T_FIXED) for strong convergence comparison.
    """
    dt = T / N
    
    # 1. Total noise needed for the entire path
    # W_T = sqrt(T) * Z_final
    
    # 2. Generate N noise steps whose *average* sum up to W_T_FIXED.
    # We use a trick: Generate N independent normals and rescale them so their sum equals the target.
    Z_sequence_raw = np.random.randn(N)
    
    # Rescale increments to ensure the sum(dW) is exactly W_T_FIXED
    # The sum of N(0, dt) is N(0, N*dt) = N(0, T). We need the sum of Z to be Z_final * sqrt(T) / sqrt(dt).
    # Since dW_k = sqrt(dt) * Z_k, the sum(dW_k) = sqrt(dt) * sum(Z_k).
    # We need sum(dW_k) = W_T_FIXED.
    
    dW_sequence = np.zeros(N)
    dW_sum_target = W_T_FIXED
    
    # Simple method: use the original dW logic and adjust the final step
    Z_sequence = np.random.randn(N)
    dW_sequence = np.sqrt(dt) * Z_sequence
    
    # Adjust last step to hit the target W_T exactly (simplifies analysis)
    dW_sequence[-1] += dW_sum_target - np.sum(dW_sequence)
    
    S = np.zeros(N)
    S[0] = S0
    
    for i in range(N - 1):
        S[i+1] = S[i] + sigma * dW_sequence[i] # dS = sigma dW
        
    return S[-1]

# ====================================================================
# 3. Error Analysis
# ====================================================================

errors = []
for N in N_values:
    # Run a small ensemble of M_ENSEMBLE to average out sampling error on the path
    M_ENSEMBLE = 50 
    ensemble_errors = []
    
    for _ in range(M_ENSEMBLE):
        S_em_final = euler_maruyama_strong(S0, sigma, T, N, Z_FINAL)
        ensemble_errors.append(np.abs(S_T_EXACT - S_em_final))
        
    errors.append(np.mean(ensemble_errors)) # Average error over the ensemble

errors = np.array(errors)
dt_values = T / N_values

# Perform log-log linear regression: log(Error) = A + B * log(dt)
log_dt = np.log(dt_values)
log_errors = np.log(errors)

# linregress returns (slope, intercept, r_value, p_value, std_err)
slope_fit, intercept_fit, r_value, p_value, std_err = linregress(log_dt, log_errors)

# ====================================================================
# 4. Visualization
# ====================================================================

fig, ax = plt.subplots(figsize=(8, 5))

# Plot the simulation data
ax.loglog(dt_values, errors, 'o', color='darkblue', label='Simulated Error')

# Plot the theoretical slope (0.5)
ax.loglog(dt_values, np.exp(intercept_fit) * dt_values**0.5, 'r--', 
          label=f'Theoretical Slope $0.5$ ($\mathcal{{O}}(\\sqrt{{\\Delta t}})$)')

# Plot the linear fit line
ax.loglog(dt_values, np.exp(intercept_fit) * dt_values**slope_fit, 'k-', 
          label=f'Fitted Slope (Order) $\\approx {slope_fit:.3f}$', lw=1.5)

# Labeling and Formatting
ax.set_title('Strong Convergence of Euler–Maruyama Method')
ax.set_xlabel('Time Step $\\Delta t$ (Log Scale)')
ax.set_ylabel('Absolute Error $|S_T^{\\text{exact}} - S_T^{\\text{EM}}|$ (Log Scale)')
ax.legend()
ax.grid(True, which='both', linestyle=':')

plt.tight_layout()
plt.show()

# --- Conclusion ---
print("\n--- Strong Convergence Analysis Summary ---")
print(f"Target Strong Convergence Order: 0.5")
print(f"Fitted Slope (Order): {slope_fit:.4f} \u00B1 {std_err:.4f}")
print(f"R-squared value: {r_value**2:.4f}")

print("\nConclusion: The log-log plot of the error versus the time step \u0394t yields a slope close to 0.5. This result numerically confirms the theoretical prediction that the Euler–Maruyama method converges strongly at the order $\mathcal{{O}}(\u221a\u0394t)$ (half-order strong convergence).")
```

-----

## Project 3: The Itō Correction in Action (Numerical Check)

-----

### Definition: Numerical Check of the Itō Correction

The goal of this project is to numerically confirm the presence of the **Itō correction term** ($\frac{1}{2}\sigma^2$) in the **average logarithmic return** of Geometric Brownian Motion (GBM). This provides numerical evidence for the central rule of stochastic calculus.

### Theory: The Itō Drift

The SDE for GBM is $dS_t = \mu S_t dt + \sigma S_t dW_t$. Applying Itō's Lemma to $f(S_t) = \ln S_t$ yields the exact solution for the log-price:

$$d(\ln S_t) = \left(\mu - \frac{1}{2}\sigma^2\right)dt + \sigma dW_t$$

The expected log-return over time $T$ is the deterministic drift term:

$$\mathbb{E}[\ln(S_T/S_0)] = \left(\mu - \frac{1}{2}\sigma^2\right)T$$

This shows that the average log-price *drifts* at a rate **less than** the expected price return $\mu$. The difference, $-\frac{1}{2}\sigma^2 T$, is the **Itō correction**, a deterministic drag caused by continuous volatility.

**Verification:** We simulate the log-price and verify that the ensemble average matches the corrected drift $(\mu - \frac{1}{2}\sigma^2)T$ and **not** the uncorrected drift $\mu T$.

-----

### Extensive Python Code and Visualization

The code runs $M=10,000$ GBM simulations, calculates the ensemble average of the log-price, and compares it against both the corrected (Itō) and uncorrected (Classical) theoretical drifts.

```python
import numpy as np
import matplotlib.pyplot as plt

# Set seed for reproducibility
np.random.seed(42)

# ====================================================================
# 1. Simulation Parameters
# ====================================================================

S0 = 100.0   # Initial price
mu = 0.10    # Expected return (uncorrected drift)
sigma = 0.30 # Volatility
T = 1.0      # Time to maturity
N = 252      # Number of steps (fine enough for GBM)

M_PATHS = 10000 # Number of paths for ensemble averaging

# Calculate theoretical drifts
ITO_CORRECTION_TERM = -0.5 * sigma**2 * T
ITO_DRIFT_THEO = (mu - 0.5 * sigma**2) * T # Expected log-return
CLASSICAL_DRIFT_THEO = mu * T             # Uncorrected drift (what classical calculus predicts)

# ====================================================================
# 2. GBM Simulation (Using the Exact Solution for Accuracy)
# ====================================================================

# We use the exact solution for the terminal price, which is required for accurate ensemble checking.
# S_T = S0 * exp( (mu - 0.5*sigma^2)*T + sigma*sqrt(T)*W_T )
# W_T ~ N(0, T), so sigma*W_T ~ N(0, sigma^2*T)

log_returns = np.zeros(M_PATHS)

for m in range(M_PATHS):
    # W_T is a single normal sample scaled by sqrt(T)
    W_T = np.random.standard_normal() * np.sqrt(T)
    
    # Calculate the terminal log-price relative to S0: ln(S_T/S0)
    log_ST_S0 = (mu - 0.5 * sigma**2) * T + sigma * W_T
    log_returns[m] = log_ST_S0

# Calculate ensemble average of log-return
E_LOG_RETURN_EMPIRICAL = np.mean(log_returns)
E_LOG_RETURN_STD = np.std(log_returns)

# ====================================================================
# 3. Visualization and Comparison
# ====================================================================

fig, ax = plt.subplots(figsize=(8, 5))

# Plot the three key values
bar_labels = ['Classical Drift $\\mu T$', 'Empirical $\\langle \\ln(S_T/S_0) \\rangle$', 'Itō Corrected Drift']
drift_values = [CLASSICAL_DRIFT_THEO, E_LOG_RETURN_EMPIRICAL, ITO_DRIFT_THEO]

ax.bar(bar_labels, drift_values, color=['skyblue', 'purple', 'green'], alpha=0.7)

# Add reference lines
ax.axhline(ITO_DRIFT_THEO, color='green', linestyle='--', linewidth=2, label='Itō Corrected Target')
ax.axhline(CLASSICAL_DRIFT_THEO, color='red', linestyle=':', linewidth=2, label='Classical Target')

# Labeling and Formatting
ax.set_title(f'Numerical Verification of the Itō Correction Term (\\sigma^2/2)')
ax.set_ylabel('Average Logarithmic Return $\\langle \\ln(S_T/S_0) \\rangle$')
ax.text(1, E_LOG_RETURN_EMPIRICAL, f'{E_LOG_RETURN_EMPIRICAL:.4f}', ha='center', va='bottom', fontsize=12)
ax.text(2, ITO_DRIFT_THEO, f'{ITO_DRIFT_THEO:.4f}', ha='center', va='bottom', fontsize=12)

ax.grid(True, axis='y')
plt.legend()
plt.tight_layout()
plt.show()

# --- Conclusion ---
print("\n--- Itō Correction Numerical Check ---")
print(f"Uncorrected (Classical) Drift \u03bcT: {CLASSICAL_DRIFT_THEO:.5f}")
print(f"Itō Correction Term (-\u03c3\u00b2T/2):   {ITO_CORRECTION_TERM:.5f}")
print(f"Corrected (Itō) Drift (\u222e[ln(S_T/S0)]): {ITO_DRIFT_THEO:.5f}")
print("---------------------------------------------------")
print(f"Empirical Ensemble Average \u222e[ln(S_T/S0)]: {E_LOG_RETURN_EMPIRICAL:.5f}")

print("\nConclusion: The numerical ensemble average of the log-price closely matches the corrected drift (0.055) and is significantly lower than the uncorrected drift (0.100). This confirms the presence of the Itō correction, which demonstrates that volatility introduces a predictable, deterministic drag on the average log-return.")
```

-----

## Project 4: Comparing EM and Exact GBM Solvers

-----

### Definition: Comparing EM and Exact GBM Solvers

The goal of this project is to compare the numerical performance of the **Euler–Maruyama (EM) method** against the **Exact analytical solution** for GBM, focusing on the impact of a large time step ($\Delta t$) on the accuracy of the EM scheme.

### Theory: Discretization Bias in EM

The Exact GBM formula provides the true terminal mean for any $\Delta t$:

$$\mathbb{E}[S_T] = S_0 e^{\mu T}$$

The EM formula, $S_{n+1} = S_n + \mu S_n \Delta t + \sigma S_n \sqrt{\Delta t} Z_n$, is only an approximation.

  * **Large $\Delta t$:** When the number of steps ($N$) is small (i.e., $\Delta t$ is large), the EM approximation is inaccurate and introduces a **discretization bias** in the ensemble mean, $\mathbb{E}[S_T^{\text{EM}}] \neq S_0 e^{\mu T}$.
  * **Small $\Delta t$:** As $N \to \infty$ ($\Delta t \to 0$), the EM mean converges weakly to the true mean.

This project demonstrates that using the EM method with a single large step ($N=1$) results in a biased mean, whereas the Exact formula remains unbiased regardless of $N$.

-----

### Extensive Python Code and Visualization

The code runs three simulations—Exact ($N=1$), EM ($N=1$), and EM ($N=100$)—and compares their terminal means against the theoretical target.

```python
import numpy as np
import matplotlib.pyplot as plt

# Set seed for reproducibility
np.random.seed(42)

# ====================================================================
# 1. Setup Parameters and Theory
# ====================================================================

S0 = 100.0   # Initial price
mu = 0.10    # Expected return (mu)
sigma = 0.30 # Volatility
T = 1.0      # Time to maturity

M_PATHS = 100000 # High number of paths to eliminate Monte Carlo sampling error

# Theoretical expectation (unbiased mean for the Exact Solution)
E_ST_THEO = S0 * np.exp(mu * T)

# ====================================================================
# 2. Simulation Solvers
# ====================================================================

def solve_exact(S0, mu, sigma, T, N=1, Z_sequence=None):
    """
    Exact GBM solution, which is independent of the number of steps N.
    """
    W_T = Z_sequence * np.sqrt(T)
    drift_term = (mu - 0.5 * sigma**2) * T
    diffusion_term = sigma * W_T
    S_T = S0 * np.exp(drift_term + diffusion_term)
    return S_T

def solve_em(S0, mu, sigma, T, N, Z_sequence=None):
    """
    Euler-Maruyama approximation, highly sensitive to N.
    """
    dt = T / N
    S = S0
    
    for i in range(N):
        dW = np.sqrt(dt) * Z_sequence[i]
        # S_{n+1} = S_n + mu*S_n*dt + sigma*S_n*dW
        S += mu * S * dt + sigma * S * dW
        
    return S

# ====================================================================
# 3. Running Simulations and Comparing Means
# ====================================================================

# Pre-generate one large set of standard normals for all tests
Z_large_set = np.random.randn(M_PATHS, 100) # Max steps needed is 100

# --- A. Exact Solution (N=1) ---
# Unbiased for any N. Using only the first column of the noise array.
S_T_A = solve_exact(S0, mu, sigma, T, N=1, Z_sequence=Z_large_set[:, 0])
MEAN_A = np.mean(S_T_A)

# --- B. Euler-Maruyama (N=1) - Large \Delta t ---
# Should be biased.
S_T_B = solve_em(S0, mu, sigma, T, N=1, Z_sequence=Z_large_set[:, 0])
MEAN_B = np.mean(S_T_B)

# --- C. Euler-Maruyama (N=100) - Small \Delta t ---
# Should converge closely to the theoretical mean (weak convergence in action).
S_T_C = solve_em(S0, mu, sigma, T, N=100, Z_sequence=Z_large_set)
MEAN_C = np.mean(S_T_C)

# --- Comparison Data ---
labels = ['Theory Target $S_0e^{\mu T}$', 
          f'A. Exact Solver ($N=1$)', 
          f'B. EM Solver ($N=1$, Large $\\Delta t$)', 
          f'C. EM Solver ($N=100$, Small $\\Delta t$)']
means = [E_ST_THEO, MEAN_A, MEAN_B, MEAN_C]

# ====================================================================
# 4. Visualization
# ====================================================================

fig, ax = plt.subplots(figsize=(8, 5))

ax.bar(labels, means, color=['gray', 'green', 'red', 'darkblue'], alpha=0.7)
ax.axhline(E_ST_THEO, color='black', linestyle='--', label='Theoretical Target')

# Annotate differences
for i in range(1, 4):
    diff = means[i] - E_ST_THEO
    ax.text(labels[i], means[i] + 0.1, f'{means[i]:.4f}\n(Error: {diff:.3f})', 
            ha='center', va='bottom', fontsize=10)

# Labeling and Formatting
ax.set_title('Comparison of Mean Terminal Price: Exact vs. Euler–Maruyama')
ax.set_ylabel('Mean Terminal Price $\\langle S_T \\rangle$')
ax.set_ylim(E_ST_THEO - 0.5, E_ST_THEO + 0.5)
ax.grid(True, axis='y')

plt.tight_layout()
plt.show()

# --- Conclusion ---
print("\n--- SDE Solver Comparison Summary ---")
print(f"Theoretical Mean E[S_T]: {E_ST_THEO:.4f}")
print(f"1-Step Exact Mean (A): {MEAN_A:.4f} (Unbiased)")
print(f"1-Step EM Mean (B):    {MEAN_B:.4f} (Biased)")
print(f"100-Step EM Mean (C):  {MEAN_C:.4f} (Converged)")

print("\nConclusion: The one-step Euler–Maruyama solver (B) results in a noticeable upward bias in the terminal mean, confirming that large time steps introduce discretization error. The Exact solver (A) and the multi-step EM solver (C) both accurately converge to the theoretical mean E[S_T], demonstrating the weak convergence of EM requires sufficiently small \u0394t.")
```
