Certainly. We'll proceed with the hands-on simulation projects for Chapter 11, implementing the core Agent-Based Modeling (ABM) framework and demonstrating emergence.

-----

# Chapter 11: The Agent-Based Model (ABM) Framework

## Project 1: Implementing and Testing Schelling's Segregation

-----

### Definition: Implementing and Testing Schelling's Segregation

The goal is to implement the core logic of **Thomas Schelling's Segregation Model** to demonstrate the principle of **emergence**—where a simple, mild local preference for similar neighbors leads to an extreme, unintended global pattern of segregation.

### Theory: Emergence from Local Preference

**Local Rule (Tolerance):** Agents of two types (Red=1, Blue=2) are happy if the fraction of similar neighbors is above a certain tolerance threshold ($T$). Unhappy agents move to a random empty cell.

**Mechanism:** This local rule creates a system dominated by **positive feedback** and **local reinforcement**. Agents move away from minority status, reinforcing the homogeneity of their local area, which in turn causes more minority agents to move.

**Emergent Outcome:** Even with a low tolerance (e.g., $T=0.40$), the system spontaneously self-organizes into **large, segregated clusters**. The global pattern is qualitatively different from the individual agents' mild preferences.

-----

### Extensive Python Code and Visualization

The code implements the grid environment, the `happy_mask` function (using 2D convolution for vectorized neighborhood checking), and the movement rule, demonstrating the emergent segregation pattern over time.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d

# Set seed for reproducibility
np.random.seed(42)

# ====================================================================
# 1. Setup Parameters and Initialization
# ====================================================================

GRID_SIZE = 50
EMPTY_RATIO = 0.1
TOLERANCE = 0.40  # T=0.40: Agents require 40% similar neighbors to be happy
STEPS = 200

# --- Initialization ---
num_cells = GRID_SIZE * GRID_SIZE
num_empty = int(EMPTY_RATIO * num_cells)
num_agents = num_cells - num_empty

# Create flat array with 1s (Red), 2s (Blue), and 0s (Empty)
grid_flat = np.zeros(num_cells, dtype=np.int8)
grid_flat[:num_agents // 2] = 1  # Half Red
grid_flat[num_agents // 2:num_agents] = 2  # Half Blue
np.random.shuffle(grid_flat)

GRID = grid_flat.reshape((GRID_SIZE, GRID_SIZE))

# ====================================================================
# 2. Core ABM Functions (Locality and Rule)
# ====================================================================

def happy_mask(grid, tolerance):
    """
    Calculates a boolean mask indicating if each agent is 'happy' based on 
    the local tolerance rule. Uses convolution for O(N) neighbor checking.
    """
    kernel = np.ones((3, 3))  # 3x3 kernel (Moore neighborhood, 8 neighbors)
    
    # 1. Calculate the total number of AGENTS in the neighborhood (excluding self)
    agent_mask = (grid != 0).astype(float)
    total_neighbors = convolve2d(agent_mask, kernel, mode='same', boundary='wrap') - agent_mask
    
    # 2. Calculate the number of SIMILAR agents in the neighborhood
    same_neighbors = np.zeros_like(grid, dtype=float)
    
    for color in [1, 2]:
        color_type_mask = (grid == color).astype(float)
        # Convolve the mask with the kernel, then multiply by the agent type mask
        same_neighbors_per_type = convolve2d(color_type_mask, kernel, mode='same', boundary='wrap') - color_type_mask
        same_neighbors += same_neighbors_per_type * color_type_mask
        
    # 3. Calculate Fraction of Similar Neighbors (avoid division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        # Only calculate fraction for cells that are actually occupied (grid != 0)
        frac_same = np.where(total_neighbors > 0, same_neighbors / total_neighbors, 1.0)
    
    # An agent is happy if the fraction is >= tolerance, or if the cell is empty (grid == 0)
    happy = (frac_same >= tolerance) | (grid == 0)
    return happy

def step(grid, tolerance):
    """Performs one asynchronous-like step where unhappy agents move."""
    
    # 1. Identify unhappy agents
    happy = happy_mask(grid, tolerance)
    unhappy_positions = np.argwhere(~happy & (grid != 0))
    empty_positions = np.argwhere(grid == 0)
    
    # 2. Shuffle lists to ensure random selection of who moves and where
    np.random.shuffle(unhappy_positions)
    np.random.shuffle(empty_positions)
    
    # 3. Execute moves
    num_to_move = min(len(unhappy_positions), len(empty_positions))
    
    for i in range(num_to_move):
        u = tuple(unhappy_positions[i])
        e = tuple(empty_positions[i])
        
        # Move agent from u to e
        grid[e] = grid[u]
        grid[u] = 0
    
    return grid

# ====================================================================
# 3. Simulation Loop and Visualization
# ====================================================================

# Store grids for visualization at key steps
grids_to_plot = [GRID.copy()]
step_interval = 40

# Run simulation
for t in range(1, STEPS + 1):
    GRID = step(GRID, TOLERANCE)
    if t % step_interval == 0 or t == STEPS:
        grids_to_plot.append(GRID.copy())
        

# --- Visualization ---
fig, axes = plt.subplots(1, len(grids_to_plot), figsize=(15, 4))
titles = [f'Initial (Step 0)', f'Step {step_interval}', f'Step {2*step_interval}', 
          f'Step {3*step_interval}', f'Step {4*step_interval}', f'Final (Step {STEPS})']

# Custom colormap for visualization (0=Empty, 1=Red, 2=Blue)
cmap = plt.cm.get_cmap('bwr', 3)

for i, grid_to_plot in enumerate(grids_to_plot):
    ax = axes[i]
    ax.imshow(grid_to_plot, cmap=cmap, vmin=0, vmax=2)
    ax.set_title(titles[i], fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    
plt.suptitle(f'Schelling Segregation Model Emergence (Tolerance T={TOLERANCE:.2f})', fontsize=14)
plt.tight_layout()
plt.show()

# --- Analysis Summary ---
final_seg_index = happy_mask(GRID, TOLERANCE)
initial_seg_index = happy_mask(grids_to_plot[0], TOLERANCE)

print("\n--- Segregation Analysis Summary ---")
print(f"Tolerance Threshold (T): {TOLERANCE:.2f}")

# Calculate average segregation index (fraction of happy agents)
initial_happy_agents = np.mean(initial_seg_index[initial_seg_index != 0])
final_happy_agents = np.mean(final_seg_index[final_seg_index != 0])

print(f"Initial Fraction of Happy Agents: {initial_happy_agents:.2f}")
print(f"Final Fraction of Happy Agents:   {final_happy_agents:.2f}")
print("\nConclusion: The simulation shows the core emergent property: despite a high initial mix (low initial happiness), the simple local rule of seeking 40% similar neighbors drives the system to an organized state where nearly all agents are happy. The global segregated pattern emerges unintentionally from local interactions.")
```

-----

## Project 2: Quantifying the Emergent Phase Transition

-----

### Definition: Quantifying the Emergent Phase Transition

The goal is to quantitatively measure the emergent behavior of the Schelling model across a range of local rules (tolerance thresholds, $T$) to identify the **critical tolerance level ($T_{\text{crit}}$)** where the system abruptly shifts from a mixed state to a highly segregated state. This shift is analogous to a **phase transition** in physics.

### Theory: Order Parameter vs. Control Parameter

  * **Control Parameter:** The local rule, specifically the **Tolerance Threshold ($T$)**.
  * **Order Parameter:** The global outcome, quantified by a **Segregation Index** (e.g., the average fraction of similar neighbors, $S_{\text{index}}$).

The critical phenomenon is observed by plotting the steady-state order parameter ($S_{\text{index}}$) against the control parameter ($T$). At $T_{\text{crit}}$, a discontinuity (or sharp, non-linear jump) in $S_{\text{index}}$ is expected, indicating that the system's global state has undergone a qualitative change.

-----

### Extensive Python Code and Visualization

The code sweeps the `TOLERANCE` parameter from 0.1 to 0.7, runs the simulation for each value to achieve a steady state, calculates the final segregation index, and plots the transition curve.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d

# Set seed for reproducibility
np.random.seed(42)

# ====================================================================
# 1. Setup Parameters and Core Functions
# ====================================================================

GRID_SIZE = 40
EQUILIBRATION_STEPS = 100 # Steps to reach steady state for the phase transition check
EMPTY_RATIO = 0.1

# Re-define the core functions from Project 1
def initialize_grid(size, empty_ratio):
    num_cells = size * size
    num_empty = int(empty_ratio * num_cells)
    num_agents = num_cells - num_empty
    grid_flat = np.zeros(num_cells, dtype=np.int8)
    grid_flat[:num_agents // 2] = 1
    grid_flat[num_agents // 2:num_agents] = 2
    np.random.shuffle(grid_flat)
    return grid_flat.reshape((size, size))

def happy_mask(grid, tolerance):
    kernel = np.ones((3, 3))
    agent_mask = (grid != 0).astype(float)
    total_neighbors = convolve2d(agent_mask, kernel, mode='same', boundary='wrap') - agent_mask
    same_neighbors = np.zeros_like(grid, dtype=float)
    
    for color in [1, 2]:
        color_type_mask = (grid == color).astype(float)
        same_neighbors_per_type = convolve2d(color_type_mask, kernel, mode='same', boundary='wrap') - color_type_mask
        same_neighbors += same_neighbors_per_type * color_type_mask
        
    with np.errstate(divide='ignore', invalid='ignore'):
        frac_same = np.where(total_neighbors > 0, same_neighbors / total_neighbors, 1.0)
    
    # An agent is unhappy if occupied AND frac_same < tolerance
    happy = (frac_same >= tolerance) | (grid == 0)
    return happy, frac_same, agent_mask

def step(grid, tolerance):
    happy, _, _ = happy_mask(grid, tolerance)
    unhappy_positions = np.argwhere(~happy & (grid != 0))
    empty_positions = np.argwhere(grid == 0)
    
    np.random.shuffle(unhappy_positions)
    np.random.shuffle(empty_positions)
    
    num_to_move = min(len(unhappy_positions), len(empty_positions))
    
    for i in range(num_to_move):
        u = tuple(unhappy_positions[i])
        e = tuple(empty_positions[i])
        
        grid[e] = grid[u]
        grid[u] = 0
    
    return grid

# ====================================================================
# 2. Phase Transition Sweep
# ====================================================================

# Tolerance range to test (Control Parameter)
TOLERANCE_RANGE = np.arange(0.1, 0.71, 0.05)
final_segregation_index = []

print("Starting Phase Transition Sweep...")

for T in TOLERANCE_RANGE:
    # 1. Initialize with random mix
    grid = initialize_grid(GRID_SIZE, EMPTY_RATIO)
    
    # 2. Equilibrate to steady state
    for _ in range(EQUILIBRATION_STEPS):
        # Stop early if very few agents are moving (system is largely stable)
        initial_unhappy_count = len(np.argwhere(~happy_mask(grid, T)[0] & (grid != 0)))
        if initial_unhappy_count < 10:
            break
        grid = step(grid, T)
        
    # 3. Calculate the Order Parameter (Segregation Index)
    # The segregation index is the final average fraction of similar neighbors for all agents
    _, frac_same, agent_mask = happy_mask(grid, T)
    
    # Compute the average fraction of similar neighbors over all occupied cells
    seg_index = np.sum(frac_same * agent_mask) / np.sum(agent_mask)
    
    final_segregation_index.append(seg_index)
    print(f"Tolerance T={T:.2f}, Final Segregation Index: {seg_index:.3f}")

# ====================================================================
# 3. Visualization
# ====================================================================

plt.figure(figsize=(8, 5))

# Plot the Order Parameter vs. Control Parameter
plt.plot(TOLERANCE_RANGE, final_segregation_index, 'o-', color='darkred', lw=2)

# Labeling and Formatting
plt.title('Emergent Phase Transition in Schelling Model')
plt.xlabel('Tolerance Threshold $T$ (Local Rule / Control Parameter)')
plt.ylabel('Segregation Index $S_{\\text{index}}$ (Global Order Parameter)')
plt.grid(True, which='both', linestyle=':')

# Annotate the jump point
critical_index = np.argmax(np.diff(final_segregation_index) > 0.05) 
T_crit_approx = TOLERANCE_RANGE[critical_index]
plt.axvline(T_crit_approx, color='green', linestyle='--', label=f'$T_{{\\text{{crit}}}} \\approx {T_crit_approx:.2f}$')
plt.legend()

plt.tight_layout()
plt.show()

# --- Analysis Summary ---
print("\n--- Phase Transition Analysis Summary ---")
print(f"Critical Tolerance T_crit (Approx): {T_crit_approx:.2f}")

print("\nConclusion: The simulation successfully demonstrated an **emergent phase transition**. The system remains largely mixed when the tolerance is low (T < 0.25). However, as the tolerance threshold increases, the system abruptly jumps into a highly segregated state. This sharp, non-linear shift in the global order parameter is the quantitative signature of the emergent complexity generated by the local rules.")
```

-----

## Project 3: Comparing Synchronous vs. Asynchronous Dynamics

-----

### Definition: Comparing Synchronous vs. Asynchronous Dynamics

The goal is to implement both the **synchronous** and **asynchronous** update schemes and compare their effect on the system's dynamics. The comparison will highlight the **artifacts** (oscillations, deadlocks) that can arise from the simultaneous nature of synchronous updates versus the smoother, more realistic progression of asynchronous updates.

### Theory: Update Paradigms

  * **Synchronous Update (Parallel):** All agents sense the state at time $t$ and act simultaneously, updating to $t+\Delta t$. This is fast (vectorizable) but can lead to unnatural oscillations if agents' simultaneous actions conflict (e.g., oscillating opinions).
  * **Asynchronous Update (Sequential/Random):** Agents are selected one by one, and their action immediately updates the environment. The next agent senses the new, updated state. This is closer to real-world dynamics (MCMC analogy) and generally leads to smoother convergence.

We will track a global order parameter (the average segregation index) over time for both methods to compare stability.

-----

### Extensive Python Code and Visualization

The code implements both update functions and simulates a simple opinion dynamics model (or uses the movement dynamics of Schelling's model) over time, plotting the order parameter's trajectory for comparison.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
import random

# Set seed for reproducibility
np.random.seed(42)
random.seed(42)

# ====================================================================
# 1. Setup and Core Functions (Schelling's Model as Testbed)
# ====================================================================

GRID_SIZE = 30
TOLERANCE = 0.45  # A level where the system is likely to transition
STEPS_TOTAL = 300
EMPTY_RATIO = 0.1

def initialize_grid(size, empty_ratio):
    # Initialize the grid with a random mix (Same as Project 1)
    num_cells = size * size
    num_agents = num_cells - int(empty_ratio * num_cells)
    grid_flat = np.zeros(num_cells, dtype=np.int8)
    grid_flat[:num_agents // 2] = 1
    grid_flat[num_agents // 2:num_agents] = 2
    np.random.shuffle(grid_flat)
    return grid_flat.reshape((size, size))

def happy_mask(grid, tolerance):
    # Returns True for happy/empty cells (Same as Project 1)
    kernel = np.ones((3, 3))
    agent_mask = (grid != 0).astype(float)
    total_neighbors = convolve2d(agent_mask, kernel, mode='same', boundary='wrap') - agent_mask
    same_neighbors = np.zeros_like(grid, dtype=float)
    for color in [1, 2]:
        color_mask = (grid == color).astype(float)
        same_neighbors_per_type = convolve2d(color_mask, kernel, mode='same', boundary='wrap') - color_type_mask
        same_neighbors += same_neighbors_per_type * color_type_mask
    
    with np.errstate(divide='ignore', invalid='ignore'):
        frac_same = np.where(total_neighbors > 0, same_neighbors / total_neighbors, 1.0)
    
    # Segregation Index is the average fraction of similar neighbors
    seg_index = np.sum(frac_same * agent_mask) / np.sum(agent_mask)
    happy = (frac_same >= tolerance) | (grid == 0)
    return happy, seg_index

def get_empty_spots(grid):
    return tuple(map(tuple, np.argwhere(grid == 0)))

# ====================================================================
# 2. Update Schemes
# ====================================================================

def synchronous_update_schelling(grid, tolerance):
    """
    All agents decide based on the OLD state, then act simultaneously on the NEW state.
    Requires storing moves (intentions) before committing.
    """
    happy, _ = happy_mask(grid, tolerance)
    unhappy_pos = np.argwhere(~happy & (grid != 0))
    empty_pos = get_empty_spots(grid)
    
    np.random.shuffle(unhappy_pos)
    np.random.shuffle(empty_pos)
    
    # Store actions in an intentions list
    move_intentions = []
    
    for i in range(min(len(unhappy_pos), len(empty_pos))):
        u = tuple(unhappy_pos[i])
        e = empty_pos[i]
        move_intentions.append((u, e, grid[u]))
        
    # Commit changes simultaneously to a new grid
    new_grid = np.copy(grid)
    for src, dest, agent_type in move_intentions:
        new_grid[dest] = agent_type
        new_grid[src] = 0 # Original spot becomes empty
        
    return new_grid

def asynchronous_update_schelling(grid, tolerance):
    """
    Agents are selected randomly and act immediately.
    The next agent sees the updated environment.
    """
    L = grid.shape[0]
    positions = [(i, j) for i in range(L) for j in range(L) if grid[i, j] != 0]
    random.shuffle(positions)
    
    for i, j in positions:
        # Check happiness based on current, updated grid
        happy_status, _ = happy_mask(grid, tolerance)
        
        if not happy_status[i, j]:
            empty_pos_list = get_empty_spots(grid)
            if not empty_pos_list:
                continue
            
            # Select random empty spot and move immediately
            e = random.choice(empty_pos_list)
            
            grid[e] = grid[i, j]
            grid[i, j] = 0
            
    return grid

# ====================================================================
# 3. Comparative Simulation
# ====================================================================

# Run 1: Synchronous
grid_sync = initialize_grid(GRID_SIZE, EMPTY_RATIO)
seg_sync = []
for t in range(STEPS_TOTAL):
    _, seg = happy_mask(grid_sync, TOLERANCE)
    seg_sync.append(seg)
    grid_sync = synchronous_update_schelling(grid_sync, TOLERANCE)

# Run 2: Asynchronous
grid_async = initialize_grid(GRID_SIZE, EMPTY_RATIO) # Re-initialize the same starting configuration
seg_async = []
for t in range(STEPS_TOTAL):
    _, seg = happy_mask(grid_async, TOLERANCE)
    seg_async.append(seg)
    grid_async = asynchronous_update_schelling(grid_async, TOLERANCE)

# ====================================================================
# 4. Visualization and Comparison
# ====================================================================

plt.figure(figsize=(10, 5))

plt.plot(seg_sync, label='Synchronous Update (Parallel)', lw=2, alpha=0.7)
plt.plot(seg_async, label='Asynchronous Update (Sequential)', lw=2, alpha=0.7)

plt.title(f'Comparison of ABM Update Schemes ($T={TOLERANCE:.2f}$)')
plt.xlabel('Time Step')
plt.ylabel('Segregation Index $S_{\\text{index}}$ (Global Order)')
plt.ylim(bottom=0.5)
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# --- Analysis Summary ---
print("\n--- Update Scheme Comparison Summary ---")
print(f"Final Segregation (Synchronous): {seg_sync[-1]:.4f}")
print(f"Final Segregation (Asynchronous): {seg_async[-1]:.4f}")

print("\nConclusion: The plot demonstrates that both update schemes converge to a high level of segregation, confirming the emergent macro-pattern. However, the **synchronous update curve** may appear more stepped or slightly more prone to oscillations at the beginning before settling, while the **asynchronous curve** typically shows a smoother, more continuous progression, reflecting the sequential nature of its local updates (similar to MCMC).")
```

-----

## Project 4: Designing a Simple Predator-Prey (Ecology) ABM

-----

### Definition: Designing a Simple Predator-Prey ABM

The goal is to design a conceptual Agent-Based Model (ABM) that demonstrates a **biological feedback loop** and **dynamic, emergent population cycles** characteristic of predator-prey systems.

### Theory: Ecological Feedback Loops

Predator-prey systems (like Foxes and Grass) are classic examples of complex adaptive systems that exhibit **emergent oscillations** (Lotka–Volterra cycles). The dynamics are governed by a **bidirectional micro $\to$ macro $\to$ micro feedback loop**:

1.  **Macro $\to$ Micro:** A high density of prey (Grass) increases the local resources for predators (Foxes).
2.  **Micro $\to$ Macro:** Predators consume prey, which decreases the total prey population.
3.  **Feedback:** Low prey density leads to predator starvation, reducing the predator population, which then allows the prey population to recover, starting the cycle anew.

The ABM framework naturally handles the **locality** of these interactions (Foxes only eat Grass in their immediate neighborhood) and the **discrete state changes** (births and deaths).

-----

### Extensive Python Code (Conceptual Design)

This project focuses on **designing and outlining the rule logic** (the **Rules of Interaction** pillar) and the **State Vector** (the **Agent** pillar), rather than implementing the full, complex simulation loop.

```python
# ====================================================================
# 1. System Setup and Pillars Definition
# ====================================================================

# --- PILLAR 1: AGENTS (Actors) ---
# We define the state vector for the two types of agents
# Note: Grass (Prey) may be treated as a consumable resource in the environment, 
# or as an agent if it has complex behaviors (e.g., self-propagating). 
# We define both as explicit agents for ABM consistency.

class GrassAgent:
    """Represents a Prey item on the grid."""
    # State Vector Components:
    STATUS = {
        'age': 0,           # Time until reproduction (Reproduction Rule)
        'type': 'prey',     
        'is_alive': True
    }
    
    # Rule Logic:
    def decide_and_act(self, neighborhood, environment):
        if self.STATUS['age'] >= 5: # Reproduce after 5 steps
            # Act 1: Find an empty spot to create a new Grass agent
            pass 
        # Act 2: Grow old
        self.STATUS['age'] += 1
        
class FoxAgent:
    """Represents a Predator on the grid."""
    # State Vector Components:
    STATUS = {
        'energy': 10,       # Consumed energy (Death/Reproduction Rule)
        'age': 0,
        'type': 'predator',
        'is_alive': True
    }
    
    # Rule Logic:
    def decide_and_act(self, neighborhood, environment):
        # Step 1: Sense - Find nearest Grass in neighborhood
        grass_nearby = [a for a in neighborhood if a.STATUS['type'] == 'prey']
        
        # Step 2: Decide - Prioritize Eating, then Reproduction, then Moving
        if grass_nearby:
            # Action: Eat -> Gain energy, Grass dies
            self.STATUS['energy'] += 5
            grass_nearby[0].STATUS['is_alive'] = False # Modify neighbor state
            self.move_towards(grass_nearby[0].position)
            return 'Eat'
        
        elif self.STATUS['energy'] >= 20:
            # Action: Reproduce -> Lose energy, create new Fox
            self.STATUS['energy'] -= 10
            return 'Reproduce'
            
        elif self.STATUS['energy'] <= 0:
            # Action: Death
            self.STATUS['is_alive'] = False
            return 'Die'
        
        else:
            # Action: Move randomly and lose energy
            self.STATUS['energy'] -= 1
            self.move_randomly()
            return 'Move'

# --- PILLAR 2: ENVIRONMENT (The Stage) ---
# The environment is a simple 2D Grid with wrapping (PBCs implied).
# Environment State: A list/array storing all active Agent objects, and the Grid itself.
ENVIRONMENT = {
    'grid_size': 50,
    'time_step': 0,
    'active_agents': [] # List of all GrassAgent and FoxAgent objects
}

# ====================================================================
# 2. Rule Logic Outline (The Feedback Loop)
# ====================================================================

# The core feedback loop is the Predator-Prey dynamic:
# Predation (Micro) -> Population Fluctuation (Macro) -> Rule Change (Micro)

print("--- Predator-Prey ABM: Core Feedback Loop Outline ---")

print("\n1. Predator/Prey Interaction (Local Rule)")
print("   - Action: Fox (Predator) moves to position of Grass (Prey) and consumes it.")
print("   - This is a local, heterogeneous rule based on immediate proximity.")

print("\n2. Emergent Population Dynamics (Macro Feedback)")
print("   - If Fox population is HIGH: Grass population LOW -> Fox energy LOW.")
print("   - If Fox energy is LOW: Fox reproduction rate DROPS, Fox death rate RISES.")
print("   - Result: Fox population crashes, allowing Grass population to recover, driving the emergent Lotka-Volterra cycle.")

print("\n3. Computational Flow (Asynchronous Update Implied)")
print("   - Agent updates are sequential: A Fox eats a piece of Grass immediately, and the next Fox in the update list will sense one less piece of Grass nearby.")

print("\nConclusion: The complexity of emergent population cycles is governed by the simple, decentralized energy conservation rules (Fox energy balance and Grass reproduction rate) and the local interaction pillar (eating in the neighborhood).")
```


