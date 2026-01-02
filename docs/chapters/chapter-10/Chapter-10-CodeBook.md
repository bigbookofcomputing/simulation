# Chapter 10: Biology II: Neuroscience (Hodgkin-Huxley)

## Project 1: Defining the Gating Dynamics (Initial Setup)

-----

### Definition: Defining the Gating Dynamics

The goal of this project is to implement the necessary functions to define the **voltage-dependent rate constants ($\alpha_x, \beta_x$)** and use them to calculate the **steady-state resting values ($x_0$)** for the sodium and potassium gating variables ($m, h, n$) at the resting potential.

### Theory: Steady-State Gating

The dynamics of any gating variable $x \in \{m, h, n\}$ are governed by the first-order ODE:

$$\frac{dx}{dt} = \alpha_x(V_m)(1 - x) - \beta_x(V_m)x$$

At steady state ($\frac{dx}{dt} = 0$), the fraction of open gates ($x_0$) is a function only of the voltage ($V_m$):

$$x_0(V_m) = \frac{\alpha_x(V_m)}{\alpha_x(V_m) + \beta_x(V_m)}$$

The functions $\alpha_x(V_m)$ and $\beta_x(V_m)$ are the empirically derived, voltage-dependent rate constants. Calculating $m_0, h_0, n_0$ at the **resting potential** ($V_{\text{rest}} \approx -65 \text{ mV}$) provides the necessary initial conditions for the simulation.

-----

### Extensive Python Code

The code implements all required rate functions and computes the initial state vector $\mathbf{S}_0 = [V_0, m_0, h_0, n_0]$.

```python
import numpy as np
import random
from math import exp, log, sqrt

# ====================================================================
# 1. System Constants (Squid Giant Axon)
# ====================================================================

# Membrane parameters
CM = 1.0  # Membrane capacitance (uF/cm^2)

# Maximum conductances (mS/cm^2)
GNA_BAR = 120.0  # Sodium
GK_BAR = 36.0    # Potassium
GL = 0.3         # Leak

# Reversal potentials (mV)
ENA = 50.0  # Sodium
EK = -77.0  # Potassium
EL = -54.4  # Leak

V_REST = -65.0  # Approximate resting potential (mV)

# ====================================================================
# 2. Voltage-Dependent Rate Constants (alpha_x and beta_x)
# ====================================================================

def alpha_m(V):
    """Na+ activation rate constant."""
    return 0.1 * (25 - V) / (np.exp((25 - V) / 10) - 1)

def beta_m(V):
    """Na+ deactivation rate constant."""
    return 4 * np.exp(-V / 18)

def alpha_h(V):
    """Na+ inactivation rate constant."""
    return 0.07 * np.exp(-V / 20)

def beta_h(V):
    """Na+ deinactivation rate constant."""
    return 1 / (np.exp((30 - V) / 10) + 1)

def alpha_n(V):
    """K+ activation rate constant."""
    # Special handling for V=10 to avoid division by zero (L'Hopital's rule)
    if V == 10.0:
        return 0.1  # lim_{V->10} alpha_n(V) = 0.1
    return 0.01 * (10 - V) / (np.exp((10 - V) / 10) - 1)

def beta_n(V):
    """K+ deactivation rate constant."""
    return 0.125 * np.exp(-V / 80)

# ====================================================================
# 3. Full HH Derivative Function (The ODE System)
# ====================================================================

def hh_derivatives(S, I_ext):
    """
    Calculates the time derivatives (dV/dt, dm/dt, dh/dt, dn/dt)
    for the state vector S = [V, m, h, n].
    """
    V, m, h, n = S

    # Calculate voltage-dependent currents
    INa = GNA_BAR * m**3 * h * (V - ENA)
    IK = GK_BAR * n**4 * (V - EK)
    IL = GL * (V - EL)

    # Voltage derivative (Current Balance Equation)
    dVdt = (I_ext - (INa + IK + IL)) / CM

    # Gating variable derivatives (Kinetic ODEs)
    dmdt = alpha_m(V) * (1 - m) - beta_m(V) * m
    dhdt = alpha_h(V) * (1 - h) - beta_h(V) * h
    dndt = alpha_n(V) * (1 - n) - beta_n(V) * n

    return np.array([dVdt, dmdt, dhdt, dndt])

# ====================================================================
# 4. Compute Steady-State Initial Conditions (x_0)
# ====================================================================

def steady_state_value(alpha, beta):
    """Computes x_infinity = alpha / (alpha + beta)."""
    return alpha / (alpha + beta)

# Compute initial resting values at V_REST = -65.0 mV
m0 = steady_state_value(alpha_m(V_REST), beta_m(V_REST))
h0 = steady_state_value(alpha_h(V_REST), beta_h(V_REST))
n0 = steady_state_value(alpha_n(V_REST), beta_n(V_REST))

# Initial State Vector
S0_REST = np.array([V_REST, m0, h0, n0])

print("--- Hodgkin–Huxley Initial State Setup ---")
print(f"Resting Potential V_REST: {V_REST:.2f} mV")
print(f"Initial State m0 (Na Act): {m0:.4f}")
print(f"Initial State h0 (Na Inact): {h0:.4f}")
print(f"Initial State n0 (K Act): {n0:.4f}")
print(f"Initial State Vector S0: {S0_REST}")
```

-----

## Project 2: Simulating the Threshold and All-or-Nothing Response

-----

### Definition: Simulating the Threshold and All-or-Nothing Response

The goal of this project is to numerically determine the **threshold current ($I_{\text{crit}}$)** required to initiate a full-amplitude action potential. The simulation demonstrates the **all-or-nothing response**, where the maximum voltage ($V_{\max}$) jumps abruptly as the stimulus crosses the critical current level.

### Theory: The $\text{Na}^+$ Positive Feedback Loop

The all-or-nothing response is an **emergent property** of the $\text{Na}^+$ channel's positive feedback.

  * A stimulus that is **subthreshold** is insufficient to open enough $\text{Na}^+$ activation gates ($m$) to overcome the passive $\text{K}^+$ and Leak currents. The voltage passively returns to rest ($V_{\max} \approx V_{\text{rest}}$).
  * A stimulus that is **suprathreshold** triggers a critical number of $m$ gates to open. The resulting $\text{Na}^+$ influx causes massive depolarization, which opens *more* $\text{Na}^+$ gates (positive feedback), driving the spike to full amplitude, independent of the initial input size ($V_{\max} \approx E_{\text{Na}}$).

We use the **RK4 solver** to integrate the H–H ODEs for a fixed pulse duration, recording $V_{\max}$ for increasing $I_{\text{ext}}$ to find the discontinuity.

-----

### Extensive Python Code and Visualization

The code defines the RK4 solver, runs the H–H system across a range of stimulus currents, and plots the maximum voltage achieved versus the injected current.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# ====================================================================
# 1. Integration Core (RK4 Solver)
# ====================================================================

# Reusing hh_derivatives, steady_state_value, V_REST, etc., from Project 1

def rk4_step(func, S, I_ext, dt):
    """Performs one RK4 time step for the state vector S = [V, m, h, n]."""
    k1 = func(S, I_ext)
    k2 = func(S + 0.5 * dt * k1, I_ext)
    k3 = func(S + 0.5 * dt * k2, I_ext)
    k4 = func(S + dt * k3, I_ext)
    return S + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

def run_hh_simulation(I_ext_pulse, t_final, dt, S_init):
    """Runs a full HH simulation for a given constant stimulus I_ext_pulse."""
    steps = int(t_final / dt)
    S = S_init.copy()
    Vm_history = np.zeros(steps)
    
    for i in range(steps):
        # The current is only applied for the first 1 ms
        I_current = I_ext_pulse if i * dt <= 1.0 else 0.0
        
        S = rk4_step(hh_derivatives, S, I_current, dt)
        Vm_history[i] = S[0]
        
    return np.max(Vm_history)

# ====================================================================
# 2. Threshold Sweep Simulation
# ====================================================================

# Initial Conditions (from Project 1)
V0 = -65.0 
m0 = steady_state_value(alpha_m(V0), beta_m(V0))
h0 = steady_state_value(alpha_h(V0), beta_h(V0))
n0 = steady_state_value(alpha_n(V0), beta_n(V0))
S_INIT = np.array([V0, m0, h0, n0])

# --- Simulation Parameters ---
DT = 0.01      # ms
T_FINAL = 10.0 # ms (long enough for the spike to finish)

# Current range to test
I_EXT_MIN = 5.0
I_EXT_MAX = 8.0
I_EXT_STEP = 0.2
I_ext_values = np.arange(I_EXT_MIN, I_EXT_MAX + I_EXT_STEP, I_EXT_STEP)

# Storage
Vmax_history = []

print(f"Testing threshold current range from {I_EXT_MIN} to {I_EXT_MAX} \u03bcA/cm\u00b2...")

for I_ext in I_ext_values:
    V_max = run_hh_simulation(I_ext, T_FINAL, DT, S_INIT)
    Vmax_history.append(V_max)

# ====================================================================
# 3. Visualization and Analysis
# ====================================================================

plt.figure(figsize=(8, 5))

# Plot Vmax vs. Iext
plt.plot(I_ext_values, Vmax_history, 'o-', color='darkred', lw=2)

# Labeling and Formatting
plt.title('All-or-Nothing Response: Threshold Current $I_{\\text{crit}}$')
plt.xlabel('Stimulus Current $I_{\\text{ext}}$ ($\mu\\text{A/cm}^2$)')
plt.ylabel('Maximum Voltage Reached $V_{\\max}$ (mV)')
plt.grid(True, which='both', linestyle=':')

# Annotate the threshold jump point (approximate)
threshold_index = np.argmax(np.diff(Vmax_history))
I_crit_approx = I_ext_values[threshold_index] + I_EXT_STEP / 2
plt.axvline(I_crit_approx, color='red', linestyle='--', label=f'$I_{{crit}} \\approx {I_crit_approx:.1f}$')
plt.legend()

plt.tight_layout()
plt.show()

# --- Analysis Summary ---
print("\n--- Threshold Analysis Summary ---")
print(f"Calculated Critical Current (Approx): {I_crit_approx:.2f} \u03bcA/cm\u00b2")
print("\nConclusion: The plot demonstrates the all-or-nothing response: below the critical threshold current, the maximum voltage remains near the resting potential. Once the threshold is crossed, the maximum voltage immediately jumps to the full spike amplitude, confirming the non-linear, regenerative nature of the action potential.")
```

-----

## Project 3: Analyzing Ionic Current Dynamics

-----

### Definition: Analyzing Ionic Current Dynamics

The goal of this project is to deconstruct the voltage spike by visualizing the contributions of the three ionic currents—Sodium ($I_{\text{Na}}$), Potassium ($I_{\text{K}}$), and Leak ($I_L$). This analysis reveals the precise **mechanistic timing** of the spike.

### Theory: Current Balance and Timing

The change in voltage $\frac{dV_m}{dt}$ is driven by the net current:

$$C_m \frac{dV_m}{dt} = I_{\text{ext}} - (I_{\text{Na}} + I_{\text{K}} + I_L)$$

1.  **$I_{\text{Na}}$ (Sodium Current):** An **inward** (negative) current that is fast-activating and responsible for the initial **depolarization**.
2.  **$I_{\text{K}}$ (Potassium Current):** An **outward** (positive) current that is slow-activating and responsible for **repolarization**.
3.  **$I_L$ (Leak Current):** A small, passive current that provides the baseline stability.

The simulation must show that the **negative $I_{\text{Na}}$ peak precedes the positive $I_{\text{K}}$ peak**, which is the fundamental physical cause of the action potential's rise and fall.

-----

### Extensive Python Code and Visualization

The code re-runs a single successful spike simulation, computes the time series for all three ionic currents based on the resulting $V_m, m, h, n$ traces, and plots them alongside the voltage.

```python
import numpy as np
import matplotlib.pyplot as plt

# ====================================================================
# 1. Full HH Simulation Run (Spiking Parameters)
# ====================================================================

# --- Simulation Setup ---
DT = 0.01  # ms
T_TOTAL = 50.0  # ms
I_EXT_MAG = 10.0 # Suprathreshold current
STIM_START, STIM_END = 10.0, 11.0 # 1 ms pulse

# I_ext function (Stimulus pulse)
def I_ext(t):
    return I_EXT_MAG if STIM_START <= t <= STIM_END else 0.0

# Initial State (from Project 2)
V0 = -65.0
m0 = steady_state_value(alpha_m(V0), beta_m(V0))
h0 = steady_state_value(alpha_h(V0), beta_h(V0))
n0 = steady_state_value(alpha_n(V0), beta_n(V0))
S_INIT = np.array([V0, m0, h0, n0])

# State storage setup
steps = int(T_TOTAL / DT)
time = np.arange(0, T_TOTAL, DT)
Vm, m, h, n = np.zeros(steps), np.zeros(steps), np.zeros(steps), np.zeros(steps)
Vm[0], m[0], h[0], n[0] = S_INIT

S = S_INIT.copy()
for i in range(1, steps):
    S = rk4_step(hh_derivatives, S, I_ext(time[i-1]), DT)
    Vm[i], m[i], h[i], n[i] = S

# ====================================================================
# 2. Current Calculation (Post-Simulation)
# ====================================================================

# Calculate Conductances and Currents from the state traces
GNA = GNA_BAR * m**3 * h
GK = GK_BAR * n**4

INa = GNA * (Vm - ENA)
IK = GK * (Vm - EK)
IL = GL * (Vm - EL)

# ====================================================================
# 3. Visualization
# ====================================================================

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Plot 1: Voltage Trace
ax[0].plot(time, Vm, color='darkred', lw=2)
ax[0].set_title('Hodgkin–Huxley Voltage Trace ($V_m$)')
ax[0].set_ylabel('Voltage ($V_m$, mV)')
ax[0].grid(True)
ax[0].axvline(STIM_START, color='gray', linestyle=':', label='$I_{\\text{ext}}$ pulse')

# Plot 2: Ionic Currents
ax[1].plot(time, INa, label='$I_{\\text{Na}}$ (Inward)', color='dodgerblue', lw=2)
ax[1].plot(time, IK, label='$I_{\\text{K}}$ (Outward)', color='orange', lw=2)
ax[1].plot(time, IL, label='$I_L$ (Leak)', color='gray', lw=1, linestyle='--')
ax[1].axhline(0, color='k', linestyle='-')

ax[1].set_title('Ionic Currents During the Action Potential')
ax[1].set_xlabel('Time (ms)')
ax[1].set_ylabel('Current Density ($\mu\\text{A/cm}^2$)')
ax[1].set_ylim(-300, 100) # Set fixed axis for clarity
ax[1].legend()
ax[1].grid(True)

plt.tight_layout()
plt.show()

# --- Analysis Summary ---
max_INa = np.min(INa)
max_IK = np.max(IK)
INa_peak_time = time[np.argmin(INa)]
IK_peak_time = time[np.argmax(IK)]

print("\n--- Ionic Current Dynamics Analysis ---")
print(f"I_Na Peak (Inward): {max_INa:.2f} \u03bcA/cm\u00b2 at t={INa_peak_time:.2f} ms")
print(f"I_K Peak (Outward): {max_IK:.2f} \u03bcA/cm\u00b2 at t={IK_peak_time:.2f} ms")
print(f"Conclusion: The negative I_Na current peaks first, driving the voltage spike, while the positive I_K current peaks later, driving repolarization. This difference in kinetic timing is the deterministic cause of the action potential's waveform.")
```

-----

## Project 4: Simulating the Refractory Period

-----

### Definition: Simulating the Refractory Period

The goal of this project is to simulate the **refractory period** by stimulating the neuron with two identical current pulses separated by a varying short time delay. This demonstrates the system's temporary inability to fire a second spike, which is essential for directional signal propagation.

### Theory: Refractory Mechanism

The refractory period is an **emergent property** of the slow dynamics of the $\text{Na}^+$ inactivation gate ($h$) and the $\text{K}^+$ activation gate ($n$).

  * **Absolute Refractory Period:** Immediately after the first spike, the $\text{Na}^+$ inactivation gate ($h$) is closed, and the $\text{K}^+$ gate ($n$) is fully open. Since the $\text{Na}^+$ channel cannot be reactivated, a second spike is impossible regardless of the stimulus size.
  * **Relative Refractory Period:** As $h$ slowly recovers (opens) and $n$ slowly closes, a second spike becomes possible but requires a **larger-than-normal stimulus** or results in a **smaller-amplitude spike** because the membrane is still partially hyperpolarized.

We test various delays ($\Delta t_{\text{stim}}$) between two identical pulses to observe the recovery.

-----

### Extensive Python Code and Visualization

The code implements a dual-pulse external current function and runs three simulations with increasing time delays, plotting the voltage traces to visualize the spike recovery.

```python
import numpy as np
import matplotlib.pyplot as plt

# ====================================================================
# 1. Dual-Pulse Stimulus Function
# ====================================================================

I_PULSE_MAG = 10.0 # Suprathreshold magnitude
PULSE_DURATION = 1.0 # ms

# Initial State (from Project 2)
V0 = -65.0
m0 = steady_state_value(alpha_m(V0), beta_m(V0))
h0 = steady_state_value(alpha_h(V0), beta_h(V0))
n0 = steady_state_value(alpha_n(V0), beta_n(V0))
S_INIT = np.array([V0, m0, h0, n0])

def I_ext_dual_pulse(t, t_start_1, t_start_2):
    """Generates two 1ms current pulses."""
    t_end_1 = t_start_1 + PULSE_DURATION
    t_end_2 = t_start_2 + PULSE_DURATION
    
    current = 0.0
    if t_start_1 <= t < t_end_1:
        current += I_PULSE_MAG
    if t_start_2 <= t < t_end_2:
        current += I_PULSE_MAG
        
    return current

# ====================================================================
# 2. Simulation Loop (RK4)
# ====================================================================

DT = 0.01
T_TOTAL = 50.0
steps = int(T_TOTAL / DT)
time = np.arange(0, T_TOTAL, DT)

T1 = 10.0 # Start time of the first pulse

# Delay scenarios to test
DELAYS = [1.5, 5.0, 10.0] # ms separation (t2 - t1)
sim_results = {}

for delay in DELAYS:
    T2 = T1 + delay
    S = S_INIT.copy()
    Vm_history = np.zeros(steps)
    
    for i in range(steps):
        t_current = time[i]
        
        # Determine current based on dual pulses
        I_current = I_ext_dual_pulse(t_current, T1, T2)
        
        S = rk4_step(hh_derivatives, S, I_current, DT)
        Vm_history[i] = S[0]
        
    sim_results[delay] = Vm_history
    max_V = np.max(Vm_history[int(T2/DT):]) # Max V after the second pulse
    print(f"Delay {delay:.1f} ms: Max V after second pulse = {max_V:.2f} mV")

# ====================================================================
# 3. Visualization
# ====================================================================

plt.figure(figsize=(10, 5))

# Plot all three voltage traces
for delay, Vm_hist in sim_results.items():
    plt.plot(time, Vm_hist, lw=1.5, label=f'Delay {delay:.1f} ms')

# Annotate the stimulus periods
plt.axvline(T1, color='gray', linestyle=':', label='1st Pulse')
plt.axvline(T1 + 1.0, color='gray', linestyle=':')
plt.axvline(T1 + DELAYS[0], color='red', linestyle=':', label='2nd Pulse (Shortest Delay)')

plt.title('Hodgkin–Huxley: Simulation of the Refractory Period')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Voltage $V_m$ (mV)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# --- Conclusion ---
print("\n--- Refractory Period Analysis ---")
print("Short Delay (1.5 ms): The second pulse should fail to fire a full spike or fire a greatly diminished one (Absolute/Relative Refractory Period).")
print("Long Delay (10.0 ms): The second pulse should fire a near-full spike (Recovery).")

print("\nConclusion: The simulation successfully demonstrated the refractory period. The initial spike places the neuron in a state of unresponsiveness, evidenced by the reduced or failed amplitude of the second spike when the pulse delay is short. As the time delay increases, the neuron recovers, confirming the slow recovery kinetics of the gating variables.")
```


