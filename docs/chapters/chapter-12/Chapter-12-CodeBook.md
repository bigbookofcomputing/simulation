Certainly. We will now implement the simulation projects for Chapter 12, focusing on Agent-Based Market Models (ABMs).

-----

# Chapter 12: Finance IV: Agent-Based Market Models

## Project 1: Modeling the Ising Market Hamiltonian

-----

### Definition: Modeling the Ising Market Hamiltonian

The goal is to implement the core local interaction rules of the **Ising Market Model**. This project demonstrates how the familiar **Metropolis update** from statistical physics (Chapter 2) can be reinterpreted as a **trader's decision** process influenced by **herding ($J$)** and **uncertainty ($T$)**.

### Theory: Market Magnetization

The system is defined by the **Ising Hamiltonian**, where $s_i$ is a trader's action ($\pm 1$ for Buy/Sell), $J$ is the herding strength, and $H$ is the external news bias:

$$E(\mathbf{s}) = -J \sum_{\langle i, j \rangle} s_i s_j - H \sum_i s_i$$

**Trader Dynamics (Metropolis):** A trader flips their action ($s_i \to -s_i$) based on the energy change $\Delta E$ and the market temperature $T$ (uncertainty).

**Net Order Flow (Magnetization):** The aggregate sentiment, or **Net Order Flow ($M$)**, is the average spin:

$$M = \frac{1}{N} \sum_i s_i$$

**Hypothesis:** At **low $T$** (low uncertainty, high $\beta$), herding dominates, leading to strong consensus and high $|M|$. At **high $T$** (high uncertainty, low $\beta$), random decisions dominate, leading to a balanced market with $M \approx 0$.

-----

### Extensive Python Code and Visualization

The code implements the Ising Market simulation for $T_{\text{low}}$ (ordered/consensus market) and $T_{\text{high}}$ (chaotic/random market) and compares the resulting magnetization (Net Order Flow).

```python
import numpy as np
import matplotlib.pyplot as plt
import random

# Set seed for reproducibility
np.random.seed(42)
random.seed(42)

# ====================================================================
# 1. Setup Parameters and Core Metropolis Functions (Ising)
# ====================================================================

# --- Market Parameters ---
N = 30                     # Grid size (N x N traders)
J = 1.0                    # Herding Strength (Coupling Constant)
H = 0.0                    # External News Bias (Field)
MCS_RUN = 5000             # Monte Carlo Sweeps
EQUILIBRATION_MCS = 500

# Critical Temperature (T_c approx 2.269) separates ordered from disordered
T_CRITICAL = 2.269

# Temperature Scenarios (T = 1/beta)
T_LOW = 1.0                # T < T_c: Ordered/Consensus Market
T_HIGH = 5.0               # T > T_c: Disordered/Chaotic Market

# --- Ising Core Functions ---
def create_lattice(N, initial_state='random'):
    """Initializes the market with random Buy/Sell decisions (+1 or -1)."""
    return np.random.choice([-1, 1], size=(N, N), dtype=np.int8)

def get_local_field(i, j, spins, N=N, J=J, H=H):
    """Calculates the local influence (field) on trader (i, j) from neighbors and news."""
    # Periodic boundary conditions (PBCs) are essential
    up = spins[(i - 1) % N, j]
    down = spins[(i + 1) % N, j]
    left = spins[i, (j - 1) % N]
    right = spins[i, (j + 1) % N]
    
    # h_local = J * sum(neighbors) + H
    return J * (up + down + left + right) + H

def metropolis_update(spins, T, J=J, H=H):
    """
    Performs one full Monte Carlo Sweep (MCS) for the market dynamics.
    Trader (i, j) flips action based on the Boltzmann probability.
    """
    N = spins.shape[0]
    total_spins = N * N
    
    beta = 1.0 / T # Uncertainty parameter
    
    for _ in range(total_spins):
        # 1. Select a random trader
        i, j = random.randrange(N), random.randrange(N)
        
        # 2. Calculate the energy change for flipping the action
        h = get_local_field(i, j, spins, N, J, H)
        # dE = 2 * current_spin * h
        dE = 2 * spins[i, j] * h 
        
        # 3. Acceptance Rule (Metropolis)
        if dE < 0 or random.random() < np.exp(-dE * beta):
            spins[i, j] *= -1 # Trader flips action (Buy <-> Sell)

def calculate_magnetization(spins):
    """Calculates Net Order Flow (Magnetization) M."""
    return np.mean(spins)

# ====================================================================
# 2. Simulation and Magnetization Comparison
# ====================================================================

def run_market_simulation(T):
    spins = create_lattice(N, initial_state='random')
    M_series = []
    
    # 1. Equilibration (Thermalization)
    for _ in range(EQUILIBRATION_MCS):
        metropolis_update(spins, T)
    
    # 2. Measurement
    for _ in range(MCS_RUN):
        metropolis_update(spins, T)
        M_series.append(calculate_magnetization(spins))
        
    return np.array(M_series), spins

# --- Run Scenarios ---
M_low_T, spins_low_T = run_market_simulation(T_LOW)
M_high_T, spins_high_T = run_market_simulation(T_HIGH)

# Calculate final ensemble averages
M_avg_low = np.mean(np.abs(M_low_T)) # Use absolute M for a measure of consensus magnitude
M_avg_high = np.mean(np.abs(M_high_T))

# ====================================================================
# 3. Visualization and Analysis
# ====================================================================

fig, ax = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: Low T Magnetization (Consensus)
ax[0].plot(M_low_T, lw=1.5, color='darkred')
ax[0].axhline(0, color='gray', linestyle='--')
ax[0].set_title(f'T={T_LOW:.1f} (Low Uncertainty): Consensus')
ax[0].set_xlabel('MCS')
ax[0].set_ylabel('Net Order Flow ($M$)')
ax[0].set_ylim(-1.1, 1.1)
ax[0].text(0.05, 0.9, f'$\langle |M| \\rangle = {M_avg_low:.2f}$', transform=ax[0].transAxes)
ax[0].grid(True)

# Plot 2: High T Magnetization (Chaos)
ax[1].plot(M_high_T, lw=1.5, color='darkblue')
ax[1].axhline(0, color='gray', linestyle='--')
ax[1].set_title(f'T={T_HIGH:.1f} (High Uncertainty): Chaos')
ax[1].set_xlabel('MCS')
ax[1].set_ylabel('Net Order Flow ($M$)')
ax[1].set_ylim(-1.1, 1.1)
ax[1].text(0.05, 0.9, f'$\langle |M| \\rangle = {M_avg_high:.2f}$', transform=ax[1].transAxes)
ax[1].grid(True)

# Plot 3: Comparison
ax[2].bar(['Low T (Consensus)', 'High T (Chaos)'], [M_avg_low, M_avg_high], color=['darkred', 'darkblue'])
ax[2].set_title('Net Order Flow Magnitude Comparison $\\langle |M| \\rangle$')
ax[2].set_ylabel('Average Magnetization')
ax[2].grid(True, axis='y')

plt.tight_layout()
plt.show()

# --- Analysis Summary ---
print("\n--- Ising Market Model Analysis Summary ---")
print(f"Average Net Order Flow (Low T={T_LOW}): {M_avg_low:.4f} (High Consensus)")
print(f"Average Net Order Flow (High T={T_HIGH}): {M_avg_high:.4f} (Randomized)")

print("\nConclusion: The simulation confirms the physics analogy. At low temperature (low uncertainty), herding (J) dominates, locking the market into a strong consensus (high net order flow). At high temperature (high uncertainty), random fluctuations prevent collective alignment, resulting in a randomized market with M \u2248 0.")
```

-----

## Project 2: Simulating Price Dynamics from Net Order Flow

-----

### Definition: Simulating Price Dynamics from Net Order Flow

The goal of this project is to simulate a rudimentary price path by treating the Ising model's **Magnetization ($M$)** as the **Net Order Flow** and using it to update the price in a simple linear model: $\Delta P_t \propto M_t$. This demonstrates how collective sentiment drives the emergent random walk of prices.

### Theory: Price as a Cumulative Random Walk

In this minimal Agent-Based Market, the price evolution is the result of accumulated order imbalance:

$$P_{t+1} = P_t + \alpha M_t + \epsilon_t$$

Where:

  * $M_t$ is the Net Order Flow (magnetization) at time $t$ (derived from the Ising dynamics).
  * $\alpha$ is a price impact coefficient.
  * $\epsilon_t$ is minimal background noise.

We run the Ising simulation near $T_c$ to capture maximum fluctuation, and the resulting price path $P(t)$ is expected to exhibit **random walk-like features** (similar to Brownian motion), driven by the spontaneous, self-organized fluctuations of the traders' collective sentiment.

-----

### Extensive Python Code and Visualization

The code reuses the Ising core, runs the simulation near the critical temperature ($T=2.5$), and plots the resulting price path generated solely by the internal dynamics of the market sentiment.

```python
import numpy as np
import matplotlib.pyplot as plt
import random

# Set seed for reproducibility
np.random.seed(42)
random.seed(42)

# ====================================================================
# 1. Setup Parameters and Ising Core (from Project 1)
# ====================================================================

N = 30
J = 1.0
H = 0.0
T_SIM = 2.5                 # Near critical temperature for high fluctuation
MCS_RUN = 10000             # Long run for path analysis
EQUILIBRATION_MCS = 500

P0 = 100.0                  # Initial Price
ALPHA_IMPACT = 0.5          # Price impact coefficient (\alpha)
EPSILON_NOISE = 0.01        # Minimal background noise (\epsilon)

# Re-using simplified Metropolis update from Project 1
# (local_field, metropolis_update, calculate_magnetization)

def run_market_simulation_for_price(T):
    spins = create_lattice(N, initial_state='random')
    M_series = []
    
    # 1. Equilibration
    for _ in range(EQUILIBRATION_MCS):
        metropolis_update(spins, T)
    
    # 2. Measurement (Record M_t for the price evolution)
    for _ in range(MCS_RUN):
        metropolis_update(spins, T)
        M_series.append(calculate_magnetization(spins))
        
    return np.array(M_series)

# ====================================================================
# 2. Price Path Simulation
# ====================================================================

M_t_series = run_market_simulation_for_price(T_SIM)

Price_t_series = np.zeros(MCS_RUN + 1)
Price_t_series[0] = P0
P_current = P0

# Generate a minimal background noise sequence
background_noise = np.random.normal(0, EPSILON_NOISE, MCS_RUN)

for t in range(MCS_RUN):
    M_t = M_t_series[t]
    
    # Price Update Rule: P_{t+1} = P_t + alpha * M_t + epsilon_t
    price_change = ALPHA_IMPACT * M_t + background_noise[t]
    P_current += price_change
    Price_t_series[t + 1] = P_current

# ====================================================================
# 3. Visualization and Analysis
# ====================================================================

time_steps = np.arange(MCS_RUN + 1)

fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

# Plot 1: Emergent Price Path
ax[0].plot(time_steps, Price_t_series, lw=1.5, color='darkgreen')
ax[0].set_title(f'Emergent Price Dynamics (T={T_SIM:.2f}, Driven by Net Order Flow $M_t$)')
ax[0].set_ylabel('Asset Price $P(t)$')
ax[0].grid(True)

# Plot 2: Driving Force (Net Order Flow / Magnetization)
ax[1].plot(time_steps[:-1], M_t_series, lw=1.0, color='crimson')
ax[1].axhline(0, color='gray', linestyle='--')
ax[1].set_title('Market Sentiment (Net Order Flow $M_t$)')
ax[1].set_xlabel('Time Step')
ax[1].set_ylabel('Sentiment $M_t$')
ax[1].set_ylim(-0.15, 0.15)
ax[1].grid(True)

plt.tight_layout()
plt.show()

# --- Analysis Summary ---
price_range = np.max(Price_t_series) - np.min(Price_t_series)

print("\n--- Emergent Price Path Analysis ---")
print(f"Simulation Temperature (Uncertainty): T={T_SIM:.2f} (Near Critical)")
print(f"Price Range Generated: {price_range:.2f}")

print("\nConclusion: The simulation demonstrates that the price path evolves as a random walk, with its drift and fluctuations directly reflecting the time evolution of the emergent Net Order Flow (M_t). The large swings in price are driven by periods of heightened collective alignment (high |M|) in the underlying trader sentiment.")
```
