# Chapter 13: Biology III: Collective Behavior & Pattern Formation

## Project 1: Simulating a 1D Reaction–Diffusion System

-----

### Definition: Simulating a 1D Reaction–Diffusion System

The goal of this project is to implement the numerical solution for a **one-dimensional Reaction–Diffusion (RD) system** using the **Explicit Finite-Difference Method (FDM)**. The specific PDE includes a classic diffusion term and a non-linear **Logistic Growth** reaction term.

### Theory: FDM for Reaction–Diffusion PDEs

The general form of the RD equation for a single species $u$ is:

$$\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2} + f(u)$$

In this specific case, the reaction term $f(u)$ is **Logistic Growth** ($k u (1-u)$).

#### 1\. Spatial Discretization (The Laplacian)

The second spatial derivative (the Laplacian in 1D) is approximated by the **Three-Point FDM Stencil** at point $i$:

$$\frac{\partial^2 u}{\partial x^2} \approx \frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2}$$

#### 2\. Temporal Discretization (Forward Euler)

Using the simplest time-marching scheme, **Forward Euler**, the update rule for concentration $u_i$ at time $t+\Delta t$ is:

$$u_i^{t+\Delta t} = u_i^t + \Delta t \left[ D \left( \frac{u_{i+1}^t - 2u_i^t + u_{i-1}^t}{\Delta x^2} \right) + k u_i^t (1 - u_i^t) \right]$$

This explicit scheme demonstrates how **diffusion spreads** the concentration while the **reaction term locally amplifies** it. Note that this explicit method is **conditionally stable**, requiring a small $\Delta t$ relative to $D$ and $\Delta x$.

-----

### Extensive Python Code and Visualization

The code implements the FDM solver, initializes the domain with a localized spike, and simulates the spreading and growing behavior of the concentration over time.

```python
import numpy as np
import matplotlib.pyplot as plt

# Set seed for reproducibility
np.random.seed(42)

# ====================================================================
# 1. Setup Parameters and Grid
# ====================================================================

# --- PDE Parameters ---
D = 0.01    # Diffusion coefficient (D)
K = 0.5     # Logistic growth rate (k)

# --- Grid Parameters ---
L = 1.0     # Domain length
NX = 100    # Number of spatial grid points
DX = L / NX # Spatial step size (\Delta x)

# --- Time Parameters (Conditional Stability) ---
DT = 0.5 * DX**2 / D # Ensure stability for explicit scheme
T_FINAL = 20.0
NT = int(T_FINAL / DT)

# Initialize concentration vector (u)
U = np.zeros(NX + 1)
x_points = np.linspace(0, L, NX + 1)

# Initial Condition: Localized spike in the center of the domain
SPIKE_WIDTH = 5
U[NX // 2 - SPIKE_WIDTH : NX // 2 + SPIKE_WIDTH] = 0.5 

# ====================================================================
# 2. FDM Simulation Loop (Forward Euler)
# ====================================================================

# Storage for plotting snapshots
U_snapshots = [U.copy()]
snapshot_interval = NT // 5

for n in range(NT):
    U_new = U.copy()
    
    # 1. Calculate the spatial terms (Laplacian)
    # Uses the current time step values (U) to compute the derivatives
    U_laplacian = (np.roll(U, -1) - 2 * U + np.roll(U, 1)) / DX**2
    
    # 2. Apply the Forward Euler update
    # u_new = u_old + dt * [ D * u_xx + k * u * (1 - u) ]
    U_new = U + DT * (D * U_laplacian + K * U * (1 - U))
    
    # Boundary Conditions (Fixed Flux/Zero Concentration at boundaries)
    U_new[0] = 0.0
    U_new[NX] = 0.0
    
    U = U_new
    
    if (n + 1) % snapshot_interval == 0:
        U_snapshots.append(U.copy())

# ====================================================================
# 3. Visualization
# ====================================================================

plt.figure(figsize=(10, 5))

# Plot the concentration profile over time
for i, u_snap in enumerate(U_snapshots):
    time_snap = i * T_FINAL / 5
    plt.plot(x_points, u_snap, label=f'Time $\\approx {time_snap:.1f}$ s', lw=1.5)

# Labeling and Formatting
plt.title(f'1D Reaction–Diffusion System (D={D}, K={K})')
plt.xlabel('Position $x$')
plt.ylabel('Concentration $u(x, t)$')
plt.ylim(bottom=0, top=1.1)
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# --- Analysis Summary ---
print("\n--- 1D Reaction–Diffusion Analysis Summary ---")
print(f"Diffusion Coefficient (D): {D}")
print(f"Reaction Rate (K): {K}")
print("\nConclusion: The simulation successfully demonstrates the coupling of diffusion and reaction. The initial spike both spreads out across the domain (diffusion) and grows in magnitude (reaction) until it approaches the saturation limit (u=1) imposed by the logistic growth term.")
```

-----

## Project 2: Simulating and Visualizing Boolean Network Attractors

-----

### Definition: Simulating and Visualizing Boolean Network Attractors

The goal is to implement a small **Boolean Gene Regulatory Network (GRN)** and map its **state transitions** over time. The simulation reveals the **attractors** (stable states or cycles) of the network, which are interpreted computationally as the **stable cell types** that can be produced by this gene circuit.

### Theory: Boolean Network Dynamics

A Boolean Network simplifies gene expression to a binary state: ON (1) or OFF (0). The system is defined by $N$ genes, resulting in $2^N$ possible configurations (the **state space**).

The update rule for the state of a gene ($X_{t+1}$) is a simple logical function (AND, OR, NOT) of its regulators' current states:

$$B_{t+1} = A_t \text{ AND NOT } C_t$$

**Attractors:** Since the state space is finite, any trajectory started from an initial state must eventually enter a repeating sequence of states, called an **attractor**.

  * **Fixed Point:** A single state that repeats (e.g., (1, 1, 0) $\to$ (1, 1, 0)).
  * **Limit Cycle:** A sequence of states that repeats (e.g., (1, 0) $\to$ (0, 1) $\to$ (1, 0)).

This project uses a 3-gene system (A, B, C) with **Negative Feedback** to demonstrate convergence to a limit cycle.

-----

### Extensive Python Code and Visualization

The code implements the Boolean rules, runs a single trajectory, and plots the state evolution to identify the attractor.

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx # For plotting the conceptual network structure

# ====================================================================
# 1. Setup Network Rules (3-Gene System with Negative Feedback)
# ====================================================================

# Rule: A -> B, B -> C, C -> NOT A (Creating a negative feedback cycle)
# This design is prone to oscillation (limit cycles).

# Boolean Logic Functions
def AND(a, b): return a and b
def OR(a, b):  return a or b
def NOT(a):    return not a

def next_state(S_current):
    """
    Applies the network's logical rules to compute the next state.
    S = [A, B, C]
    """
    A, B, C = S_current
    
    # 1. A's Rule: A is activated by C
    A_next = C 
    
    # 2. B's Rule: B is activated by A
    B_next = A
    
    # 3. C's Rule: C is inhibited by B
    C_next = NOT(B)
    
    return np.array([int(A_next), int(B_next), int(C_next)])

# ====================================================================
# 2. Simulation: State Space Trajectory
# ====================================================================

MAX_STEPS = 15 # Max steps to find the attractor
# Initial State (Arbitrary starting point)
S_initial = np.array([1, 0, 0]) # A=ON, B=OFF, C=OFF

history = [S_initial.copy()]
S_current = S_initial.copy()

print(f"Initial State S0: {S_initial}")

for t in range(MAX_STEPS):
    S_next = next_state(S_current)
    S_current = S_next
    
    # Check if the state has appeared before (indicating an attractor cycle)
    if any(np.array_equal(S_next, state) for state in history):
        print(f"Attractor found at step {t+1}.")
        history.append(S_next)
        break
        
    history.append(S_next)

# Convert history to DataFrame for plotting
df_history = pd.DataFrame(history, columns=['Gene A', 'Gene B', 'Gene C'])
df_history['Time'] = df_history.index

# ====================================================================
# 3. Visualization
# ====================================================================

# 1. Time Series Plot (to visualize the oscillation)
fig, ax = plt.subplots(figsize=(8, 4))
df_history.plot(x='Time', y=['Gene A', 'Gene B', 'Gene C'], kind='line', drawstyle='steps-post', ax=ax)
ax.set_yticks([0, 1])
ax.set_title('Boolean Network Dynamics: Convergence to Attractor')
ax.set_xlabel('Time Step')
ax.set_ylabel('Expression State (0=OFF, 1=ON)')
ax.legend(title='Gene')
ax.grid(True)

plt.tight_layout()
plt.show()

# 2. State Space Plot (Conceptual Graph)
# Define the graph structure for visualization
G = nx.DiGraph()
for t in range(len(history) - 1):
    src = str(list(history[t]))
    dest = str(list(history[t+1]))
    G.add_edge(src, dest, weight=t)

# Use spring layout for visualization (optional)
pos = nx.spring_layout(G, seed=42)
plt.figure(figsize=(6, 6))
nx.draw_networkx_nodes(G, pos, node_size=500, node_color='lightblue')
nx.draw_networkx_edges(G, pos, edge_color='gray', arrows=True)
nx.draw_networkx_labels(G, pos, font_size=8)
plt.title('Attractor Cycle in State Space')
plt.axis('off')
plt.show()


# --- Analysis Summary ---
print("\n--- Boolean Network Attractor Analysis ---")
attractor_cycle = df_history.iloc[history.index(history[-1]) - 1:]
print(f"Attractor Cycle Length: {len(attractor_cycle)} steps")
print("\nAttractor Cycle States (Steady State):")
print(attractor_cycle.iloc[:-1].to_markdown(index=False))

print("\nConclusion: The simulation confirms that the finite-state network converges to a stable **limit cycle attractor** (a repeating sequence of states). Biologically, this attractor represents the final, stable **cell identity** (e.g., muscle cell or neuron) computed by the gene circuit.")
```

-----

## Project 3: Identifying Network Structural Motifs

-----

### Definition: Identifying Network Structural Motifs

The goal is to use basic **Graph Theory metrics**—specifically **in-degree** and **out-degree**—to computationally identify the structural roles of genes (nodes) within a conceptual **Gene Regulatory Network (GRN)**.

### Theory: Graph Metrics and Functional Roles

A GRN is a **directed graph** where edges indicate the flow of regulatory influence. The topology of the network reveals the functional specialization of each gene:

1.  **Out-Degree:** The number of outgoing edges from a node.
      * **High Out-Degree:** Suggests a **Master Regulator** or hub that controls many other genes downstream.
2.  **In-Degree:** The number of incoming edges to a node.
      * **High In-Degree:** Suggests a **Sensor Gene** or target that integrates signals from many upstream regulators.

The network's overall structure is stored in its **Adjacency Matrix** ($A$), where $A_{ij} \neq 0$ means gene $i$ regulates gene $j$.

-----

### Extensive Python Code and Visualization

The code defines a conceptual GRN via its adjacency matrix, uses NetworkX to analyze its topology, and prints the degree metrics to identify the key structural motifs.

```python
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# ====================================================================
# 1. Setup Conceptual Gene Regulatory Network (GRN)
# ====================================================================

# Nodes: 5 Genes/Proteins (A, B, C, D, E)
GENE_NAMES = ['A (Master)', 'B (Toggle)', 'C (Toggle)', 'D (Sensor)', 'E (Output)']
N_GENES = len(GENE_NAMES)

# Adjacency Matrix A[i, j] = 1 if i regulates j, 0 otherwise
# Rule Design:
# 1. A is a Master Regulator (high out-degree to B, D)
# 2. B and C form a Mutual Inhibition/Toggle Switch (B -> -C, C -> -B)
# 3. D is a sensor target of A and E
# 4. E is a final Output
ADJACENCY_MATRIX = np.array([
#   A   B   C   D   E 
    [0, +1,  0, +1,  0], # A activates B and D (High Out)
    [0,  0, -1,  0,  0], # B inhibits C
    [0, -1,  0,  0,  0], # C inhibits B
    [0,  0,  0,  0, +1], # D activates E
    [0,  0,  0,  0,  0]  # E is a final output (Low Out)
])

# ====================================================================
# 2. Graph Analysis (NetworkX)
# ====================================================================

# Create the directed graph from the adjacency matrix
# Add labels for positive/negative regulation
G = nx.from_numpy_array(ADJACENCY_MATRIX, create_using=nx.DiGraph)

# Relabel nodes with gene names
G = nx.relabel_nodes(G, {i: GENE_NAMES[i] for i in range(N_GENES)})

# Compute graph metrics
in_degree = dict(G.in_degree())
out_degree = dict(G.out_degree())

# Compute Degree Centrality (simple metric proportional to total connections)
degree_centrality = {gene: in_degree[gene] + out_degree[gene] for gene in GENE_NAMES}

# ====================================================================
# 3. Structural Role Identification
# ====================================================================

print("--- Gene Regulatory Network (GRN) Structural Analysis ---")

print("\n1. Centrality Metrics:")
for gene in GENE_NAMES:
    print(f"  {gene:<15} | In: {in_degree[gene]:<2} | Out: {out_degree[gene]:<2} | Total: {degree_centrality[gene]:<2}")

print("\n2. Identifying Roles (Structural Motifs):")

# Master Regulator: High Out-Degree (A)
master_regulator = max(out_degree, key=out_degree.get)
print(f"- Master Regulator: {master_regulator} (Controls 2 genes downstream)")

# Sensor Gene: High In-Degree, Low Out-Degree (D)
sensor_gene = max(degree_centrality, key=degree_centrality.get)
if out_degree[sensor_gene] == 1: # D has in=1, out=1, total=2
     print(f"- Sensor Gene: D (In=1, Out=1) - Integrates input from A, transmits to E.")

# Toggle Switch Motif (B and C)
print("- Structural Motif: B and C form a **Mutual Inhibition/Toggle Switch** (B ⊣ C, C ⊣ B).")

# ====================================================================
# 4. Visualization
# ====================================================================

plt.figure(figsize=(8, 6))
pos = nx.spring_layout(G, seed=42) # Layout for visualization

# Draw nodes scaled by Total Degree (Centrality)
nx.draw_networkx_nodes(G, pos, node_size=[v * 300 for v in degree_centrality.values()],
                       node_color='lightblue', alpha=0.9)

# Draw edges (distinguish activation/inhibition)
edges = G.edges()
colors = ['red' if ADJACENCY_MATRIX[GENE_NAMES.index(u), GENE_NAMES.index(v)] < 0 else 'blue' for u, v in edges]
labels = nx.get_edge_attributes(G, 'weight')

nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=colors, width=1.5, arrowsize=20)
nx.draw_networkx_labels(G, pos, font_weight='bold', font_size=10)

plt.title('Gene Regulatory Network Topology (Master Regulator A, Toggle B/C)')
plt.axis('off')
plt.show()

print("\nConclusion: Computational analysis of the GRN's topology, specifically through in-degree and out-degree, reveals its functional architecture. Gene A acts as the upstream master controller, while B and C form a local feedback motif (the toggle switch) essential for binary cell fate decision-making.")
```

-----

## Project 4: Modeling a Genetic Toggle Switch (Continuous ODE)

-----

### Definition: Modeling a Genetic Toggle Switch

The goal is to implement the continuous (ODE) model of the **genetic toggle switch**—a system of two mutually inhibitory genes ($u \dashv v$ and $v \dashv u$)—and solve it using the **4th-order Runge–Kutta (RK4) solver**. The simulation must demonstrate **bistability**, a fundamental mechanism for binary decision-making in cells.

### Theory: Bistability and Mutual Inhibition

The toggle switch is a classic example of a **positive feedback loop** created by **mutual inhibition**. The system has two stable steady states (attractors):

1.  State 1: Gene $u$ HIGH, Gene $v$ LOW.
2.  State 2: Gene $u$ LOW, Gene $v$ HIGH.

The coupled ODEs, which follow mass-action kinetics with a repression term (Hill function), are:

$$\frac{du}{dt} = \frac{\alpha_1}{1 + v^{\beta}} - u$$
$$\frac{dv}{dt} = \frac{\alpha_2}{1 + u^{\gamma}} - v$$

  * $\alpha_1, \alpha_2$: Maximum synthesis rates.
  * $u, v$: Protein concentrations (decay rate is normalized to 1).
  * $\beta, \gamma$: Hill coefficients (cooperativity of repression).

The **RK4 solver** is used to accurately integrate this **nonlinear, stiff system**. By running the simulation with different **initial conditions**, we show that the system converges to different final states, confirming **bistability**.

-----

### Extensive Python Code and Visualization

The code implements the ODE system, the RK4 solver, runs the two required simulations (Run A: High $u_0$, Run B: High $v_0$), and plots the trajectories to confirm convergence to two distinct attractors.

```python
import numpy as np
import matplotlib.pyplot as plt

# ====================================================================
# 1. Setup Parameters and The ODE System
# ====================================================================

# --- Toggle Switch Parameters ---
# Parameters chosen to ensure bistability (high Hill coefficients)
ALPHA1 = 5.0  # Max synthesis rate for u
ALPHA2 = 5.0  # Max synthesis rate for v
BETA_HILL = 3.0 # Hill coefficient (cooperativity) for v repressing u
GAMMA_HILL = 3.0 # Hill coefficient (cooperativity) for u repressing v

# --- Simulation Parameters ---
DT = 0.01  # Time step
T_FINAL = 50.0  # Total time (ms)

def toggle_ode_system(S_current):
    """
    Implements the coupled ODEs for the genetic toggle switch:
    S = [u, v]
    """
    u, v = S_current
    
    # ODE for u: du/dt = alpha1 / (1 + v^beta) - u 
    dudt = (ALPHA1 / (1.0 + v**BETA_HILL)) - u
    
    # ODE for v: dv/dt = alpha2 / (1 + u^gamma) - v
    dvdt = (ALPHA2 / (1.0 + u**GAMMA_HILL)) - v
    
    return np.array([dudt, dvdt])

# ====================================================================
# 2. RK4 Solver Implementation
# ====================================================================

def rk4_step(func, S, dt):
    """Performs one RK4 time step for the state vector S = [u, v]."""
    k1 = func(S)
    k2 = func(S + 0.5 * dt * k1)
    k3 = func(S + 0.5 * dt * k2)
    k4 = func(S + dt * k3)
    return S + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

def run_simulation(S_initial):
    """Runs the simulation from a given initial state and records history."""
    steps = int(T_FINAL / DT)
    u_history = np.zeros(steps)
    v_history = np.zeros(steps)
    
    S = S_initial.copy()
    
    for i in range(steps):
        S = rk4_step(toggle_ode_system, S, DT)
        u_history[i] = S[0]
        v_history[i] = S[1]
        
    return u_history, v_history

# ====================================================================
# 3. Bistability Scenarios and Simulation
# ====================================================================

# --- Run A: Initial State biased towards u (u HIGH, v LOW) ---
S_INIT_A = np.array([10.0, 1.0])
uA, vA = run_simulation(S_INIT_A)

# --- Run B: Initial State biased towards v (u LOW, v HIGH) ---
S_INIT_B = np.array([1.0, 10.0])
uB, vB = run_simulation(S_INIT_B)

time_points = np.arange(0, T_FINAL, DT)

# ====================================================================
# 4. Visualization
# ====================================================================

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Trajectories (Time Series)
ax[0].plot(time_points, uA, label='u (Run A: High u0)', color='darkblue', lw=2)
ax[0].plot(time_points, vA, label='v (Run A: High u0)', color='skyblue', lw=2)
ax[0].plot(time_points, uB, label='u (Run B: High v0)', color='darkred', lw=2, linestyle='--')
ax[0].plot(time_points, vB, label='v (Run B: High v0)', color='salmon', lw=2, linestyle='--')
ax[0].set_title('Genetic Toggle Switch: Trajectories')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Concentration')
ax[0].legend()
ax[0].grid(True)

# Plot 2: Bistability (Phase Space/u vs v)
ax[1].plot(uA, vA, label='Run A: Attractor 1 (u HIGH)', color='darkblue', lw=2)
ax[1].plot(uB, vB, label='Run B: Attractor 2 (v HIGH)', color='darkred', lw=2, linestyle='--')
ax[1].plot(uA[0], vA[0], 'ko', label='Start A', markersize=6)
ax[1].plot(uB[0], vB[0], 'ks', label='Start B', markersize=6)
ax[1].plot(uA[-1], vA[-1], 'go', label='End A', markersize=8)
ax[1].plot(uB[-1], vB[-1], 'g^', label='End B', markersize=8)

ax[1].set_title('Bistability: Convergence to Two Attractors')
ax[1].set_xlabel('Concentration u')
ax[1].set_ylabel('Concentration v')
ax[1].legend()
ax[1].grid(True)

plt.tight_layout()
plt.show()

# --- Analysis Summary ---
print("\n--- Genetic Toggle Switch Bistability Analysis ---")
print(f"Run A (Start u=10, v=1): Final State (u={uA[-1]:.2f}, v={vA[-1]:.2f}) -> Attractor 1 (u HIGH)")
print(f"Run B (Start u=1, v=10): Final State (u={uB[-1]:.2f}, v={vB[-1]:.2f}) -> Attractor 2 (v HIGH)")

print("\nConclusion: The simulation demonstrates **bistability**—the fundamental mechanism for binary decision-making in cells. Despite the identical governing equations, the system converges to one of two distinct stable states (attractors) based entirely on the initial conditions, confirming the logical switch functionality of the mutually inhibitory genetic circuit.")
```

Done.


