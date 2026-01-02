## 🧠 Chapter 14: Biology IV: Computational Neuroscience (Workbook)

The goal of this chapter is to scale up from the single-neuron dynamics (Chapter 10) to **network behavior**, showing how collective computation, memory, and pattern recognition emerge from simple, coupled neural elements, using the **Hopfield Network** as the primary model.

| Section | Topic Summary |
| :--- | :--- |
| **14.1** | Chapter Opener: The Network as a Computer |
| **14.2** | The Agents: Simplification to Integrate-and-Fire |
| **14.3** | The Physics Analogy: The Hopfield Network and Memory |
| **14.4** | The Simulation: Storing and Retrieving Patterns |
| **14.5** | Chapter Summary & End of Volume II |

***

### 14.1 The Network as a Computer

> **Summary:** The complexity of the brain comes from the **collective interaction** of billions of neurons, not from a single neuron's internal dynamics. The goal of computational neuroscience at this level is to model how simple, interconnected units can perform higher-level functions like **memory and recognition**. This is achieved by drawing an analogy between neural stability and the **low-energy ground states** of physical systems like the Ising Model.

#### Section Detail

Modeling the entire brain using the complex Hodgkin-Huxley ODEs (Chapter 10) is computationally infeasible. Therefore, analysis must focus on the emergent behavior of the network, which requires simplifying the neuron to its essential signaling state. The fundamental insight is that stability in a network, whether a gene regulatory network or a memory network, is determined by its tendency to settle into **attractors**.

#### Quiz Questions

**1. Why are complex models like Hodgkin-Huxley (H-H) impractical for simulating a network of thousands of neurons?**

* **A.** H-H does not account for the action potential.
* **B.** H-H is a stochastic model.
* **C.** **H-H is too computationally expensive to simulate for large networks**. (**Correct**)
* **D.** H-H violates the principle of local interaction.

**2. The higher-level cognitive functions of the brain, such as memory and recognition, are primarily considered emergent properties of the:**

* **A.** Complexity of the neuron's membrane potential.
* **B.** **Collective interaction of billions of neurons in a network**. (**Correct**)
* **C.** The exact formula for the $K^+$ current.
* **D.** The process of transcription.

---

#### Interview-Style Question

**Question:** Explain why the concept of **emergence** is essential for understanding the brain, as opposed to a reductionist approach that focuses on individual neuron function.

**Answer Strategy:** A reductionist approach (like H-H modeling) explains the physics of the **single-neuron signal** (the spike). However, **memory, thought, and cognition** are not properties of a single neuron; they are properties of the **network**. Emergence is essential because it explains how a simple, collective pattern (e.g., a specific set of 10,000 neurons firing simultaneously) can arise from, and stabilize across, billions of individual interactions. The whole (memory) is qualitatively different from the sum of its parts (individual spikes).

---

### 14.2 The Agents: Simplification to Integrate-and-Fire

> **Summary:** To make network modeling computationally feasible, the complex H-H neuron is simplified to the **Integrate-and-Fire** model. The key simplification is the **Binary Neuron State**, where each neuron $i$ is modeled as a simple **spin** ($s_i \in \{+1, -1\}$), representing **firing** ($+1$) or **silent** ($-1$).

#### Section Detail

The simplification shifts the focus from complex **membrane dynamics** to the fundamental **network topology and synaptic weights** ($w_{ij}$). The binary spin state makes the neural network mathematically equivalent to the **Ising Model**. In the integrate-and-fire rule, the neuron fires if its total incoming current (input) exceeds a threshold.

#### Quiz Questions

**1. The primary purpose of simplifying the neuron model from Hodgkin-Huxley to the Integrate-and-Fire model for network studies is to:**

* **A.** Introduce a refractory period.
* **B.** **Shift the computational focus from complex membrane dynamics to network connectivity**. (**Correct**)
* **C.** Introduce a continuous voltage variable.
* **D.** Model the $Cl^-$ current.

**2. In the Hopfield Network's simplified model, the state $s_i = +1$ typically represents which neural event?**

* **A.** The resting potential.
* **B.** **The neuron is firing (active)**. (**Correct**)
* **C.** The $Na^+$ inactivation gate is closed.
* **D.** The neuron is inhibited.

---

#### Interview-Style Question

**Question:** The transition from a continuous $\text{H-H}$ voltage variable ($V_m$) to a discrete binary state ($s_i = \pm 1$) is a massive abstraction. What crucial aspect of neural signaling does the binary simplification still successfully capture?

**Answer Strategy:** The binary state successfully captures the **all-or-nothing nature** of the action potential (Chapter 10). A neuron either fires at its full amplitude ($+1$) or it doesn't ($-1$). Since information in the brain is often encoded in the **rate and pattern of firing** rather than the precise voltage of a single spike, the binary abstraction retains the essential functional output needed for network-level computation.

---

### 14.3 The Physics Analogy: The Hopfield Network and Memory

> **Summary:** The **Hopfield Network** is a fully connected neural network where the network state $\mathbf{s}$ is related to an **Energy Function** $E(\mathbf{s}) = - \frac{1}{2} \sum_{i, j} w_{ij} s_i s_j$. This function is mathematically identical to the interaction term of the **Ising Hamiltonian**. **Memories** are encoded by the **synaptic weights ($w_{ij}$)** using the **Hebbian learning rule** and are stored as **stable, low-energy minima (attractors)** in the network's energy landscape.

#### Section Detail

The retrieval process in the Hopfield network is equivalent to the system undergoing a **gradient descent** (or relaxation) in the energy landscape. When the network is given a corrupted input (a "cue"), it evolves until it reaches the nearest low-energy minimum, which corresponds to the **closest stored memory**. The memory is the pattern of activation that stabilizes the system.

#### Quiz Questions

**1. In the Hopfield Network's Energy Function, $E(\mathbf{s}) = - \frac{1}{2} \sum w_{ij} s_i s_j$, what variable is analogous to the $\frac{1}{k_B T}$ term in the Ising Model?**

* **A.** The synaptic weight $w_{ij}$.
* **B.** The neuron state $s_i$.
* **C.** **The temperature (or stochastic noise) of the neural system**. (**Correct**)
* **D.** The number of stored patterns $M$.

**2. In the Hopfield Network, a stored memory corresponds to which feature in the network's energy landscape?**

* **A.** The global maximum of the energy function.
* **B.** A random walk.
* **C.** **A stable, low-energy minimum (attractor)**. (**Correct**)
* **D.** A high-energy metastable state.

---

#### Interview-Style Question

**Question:** Explain the concept of **associative memory** using the terms *attractor*, *energy landscape*, and *Hebbian learning*.

**Answer Strategy:**
1.  **Encoding (Hebbian Learning):** Patterns (memories) are encoded into the network by setting the **synaptic weights ($w_{ij}$) using the Hebbian rule**, which essentially determines the shape of the **energy landscape**.
2.  **Storage (Attractors):** This process makes the encoded patterns correspond to specific **low-energy minima (attractors)** in the energy landscape.
3.  **Retrieval (Relaxation):** When the network is given a **partial or corrupted input** (a cue), the system undergoes **relaxation** (gradient descent) and evolves from the current high-energy state until it falls into the nearest low-energy attractor, thereby **completing the pattern** and recalling the associated memory.

---

### 14.4 The Simulation: Storing and Retrieving Patterns

> **Summary:** The computational process involves two phases: **Encoding**, where the synaptic weight matrix $W$ is calculated using the **Hebbian learning rule** ($w_{ij} = \frac{1}{M} \sum_{m=1}^M s^{(m)}_i s^{(m)}_j$), and **Retrieval**, where the network is initialized with a corrupted pattern and evolves iteratively using an **asynchronous update** until it reaches a stable state (attractor).

#### Section Detail

The asynchronous update, similar to MCMC in Chapter 2, involves selecting one random neuron $i$ at a time, calculating its weighted input ($h_i$), and updating its state ($s_i^{\text{new}} = \pm 1$) based on a threshold. This process ensures the system is always moving toward a lower energy state ($ \Delta E \le 0$), confirming the memory retrieval mechanism is a relaxation process.

#### Quiz Questions

**1. The **Hebbian learning rule** is the mechanism used in the Hopfield Network to calculate the synaptic weights $W$ during the encoding phase. This rule primarily relies on:**

* **A.** The external magnetic field $H$.
* **B.** The nearest-neighbor interaction $J$.
* **C.** **The outer product (correlation) of the patterns being stored**. (**Correct**)
* **D.** The $K^+$ current dynamics.

**2. The dynamics of memory retrieval in the Hopfield Network is simulated by an iterative, asynchronous update that ensures the network is always moving toward a state where:**

* **A.** The total voltage is maximized.
* **B.** **The network energy $E(\mathbf{s})$ is minimized**. (**Correct**)
* **C.** The temperature is maximized.
* **D.** All neurons are firing ($s_i = +1$).

---

#### Interview-Style Question

**Question:** In the Hopfield Network, the asynchronous update rule is deterministic (no acceptance probability, unlike Metropolis). Why is it guaranteed that this deterministic process will eventually stabilize (stop changing) and settle into a minimum, rather than oscillating forever?

**Answer Strategy:** The update rule is constructed to be a pure **gradient descent** process. By definition, the rule only permits moves that either **lower the network energy ($\Delta E < 0$) or leave it unchanged ($\Delta E = 0$)**. Since the energy function $E(\mathbf{s})$ is bounded from below (it has a minimum possible value), the system must eventually run out of energy-lowering moves and therefore stabilize in a state where $\Delta E = 0$, which is a local minimum (the attractor).

---

## 💡 Hands-On Simulation Projects (Chapter Conclusion) 🛠️

These projects require implementing the core Hopfield Network mechanics and testing its memory properties, drawing analogies to the Ising Model.

### Project 1: Encoding and Analyzing the Weight Matrix ($W$)

* **Goal:** Implement the Hebbian learning rule and analyze the structure of the resulting weight matrix.
* **Setup:** Define three simple, orthogonal (low-overlap) binary patterns $\mathbf{s}^{(1)}, \mathbf{s}^{(2)}, \mathbf{s}^{(3)}$ on a small $N=10$ neuron network (e.g., $\mathbf{s}^{(1)} = [1, 1, \dots, 1]$).
* **Steps:**
    1.  Implement the Hebbian learning formula $w_{ij} = \frac{1}{M} \sum_{m=1}^M s^{(m)}_i s^{(m)}_j$ to calculate the $W$ matrix.
    2.  Plot the resulting $W$ matrix as a heatmap.
* ***Goal***: Show that the diagonal elements ($w_{ii}$) are positive (self-excitatory), and that the off-diagonal elements reflect the correlation between the stored patterns.

### Project 2: Simulating Pattern Retrieval and Error Correction

* **Goal:** Demonstrate the network's ability to perform **associative recall** (memory retrieval) from a corrupted cue.
* **Setup:** Use the $W$ matrix from Project 1. Define a noisy input pattern $\mathbf{s}_{\text{input}}$ by randomly flipping $20\%$ of the bits of one stored pattern $\mathbf{s}^{(1)}$.
* **Steps:**
    1.  Implement the asynchronous retrieval loop: select a random neuron $i$, calculate $h_i = \sum w_{ij} s_j$, and update $s_i^{\text{new}} = \text{sgn}(h_i)$.
    2.  Run the simulation for 100 steps.
    3.  Measure the **overlap** of the network state $\mathbf{s}(t)$ with the original, uncorrupted pattern $\mathbf{s}^{(1)}$ at each step.
* ***Goal***: Show that the overlap quickly increases from $80\%$ to $100\%$, confirming that the network relaxes toward and retrieves the full, correct stored pattern.

### Project 3: Visualizing the Energy Landscape and Relaxation

* **Goal:** Show that the memory retrieval process is a genuine **gradient descent** (relaxation) in the energy landscape.
* **Setup:** Use the same network and retrieval simulation from Project 2.
* **Steps:**
    1.  Write a function to calculate the network's total **Energy** $E(\mathbf{s})$ at any given state $\mathbf{s}$.
    2.  During the retrieval simulation, calculate and record $E(t)$ at each time step.
* ***Goal***: Plot $E(t)$ versus time. The plot should be monotonically decreasing or constant ($\Delta E \le 0$), confirming that the asynchronous update rule always drives the system downhill toward a minimum.

### Project 4: Testing Network Capacity (The Fidelity Analogy)

* **Goal:** Demonstrate the limitation of the Hopfield Network: too many stored memories cause interference (analogous to the fidelity limits in a recording system).
* **Setup:** Use a larger $N=100$ neuron network.
* **Steps:**
    1.  Run the encoding (Hebbian rule) for a small number of patterns (e.g., $M=5$). Test retrieval fidelity (Project 2).
    2.  Run the encoding again for a large number of patterns (e.g., $M=50$). Test retrieval fidelity.
* ***Goal***: Show that retrieval fidelity remains high for $M=5$, but **fails catastrophically** for $M=50$ (the network retrieves a mixture or a spurious state). This demonstrates the hard limit on memory capacity imposed by the physics of the energy landscape.

## 🧠 Chapter 14: Biology IV: Computational Neuroscience (Workbook)

The goal of this chapter is to scale up from the single-neuron dynamics (Chapter 10) to **network behavior**, showing how collective computation, memory, and pattern recognition emerge from simple, coupled neural elements, using the **Hopfield Network** as the primary model.

| Section | Topic Summary |
| :--- | :--- |
| **14.1** | Chapter Opener: The Network as a Computer |
| **14.2** | The Agents: Simplification to Integrate-and-Fire |
| **14.3** | The Physics Analogy: The Hopfield Network and Memory |
| **14.4** | The Simulation: Storing and Retrieving Patterns |
| **14.5** | Chapter Summary & End of Volume II |

***

### 14.1 The Network as a Computer

> **Summary:** The brain's intelligence is **distributed** across a complex network of neurons, not held by a single cell. The network is modeled as an **associative memory system** where memories correspond to **stable attractor states** in a high-dimensional energy landscape. This framework unifies physics and cognition.

#### Section Detail

The **Hopfield Network** is the model that formalizes this connection, realizing that the network's activity evolves to minimize a scalar **energy function**, just like a physical system. This dynamics-driven computation is a self-organizing process where **local physics yields global intelligence**.

#### Quiz Questions

**1. The primary breakthrough of the Hopfield Network model was realizing that:**

* **A.** Neuron firing is continuous.
* **B.** **The network's activity evolves to minimize a scalar energy function**. (**Correct**)
* **C.** Memories are stored in a central location.
* **D.** Only inhibitory connections are needed.

**2. In the Hopfield model, what does the network's system dynamics eventually settle into, which is interpreted as a "stored memory"?**

* **A.** A magnetic field.
* **B.** **A stable, low-energy attractor**. (**Correct**)
* **C.** The exact initial input cue.
* **D.** A Boltzmann distribution.

---

#### Interview-Style Question

**Question:** The text describes the brain's computation as closer to a **thermodynamic computer** than a **digital computer**. Explain the difference in their computational processes.

**Answer Strategy:**
* **Digital Computer:** Computes by executing sequential, step-by-step instructions based on formal logic (like a Turing Machine). The output is determined by the program.
* **Thermodynamic Computer (Hopfield/Brain):** Computes by **relaxing to equilibrium**. It starts in a high-energy state (the input cue) and evolves dynamically toward a stable, low-energy state (the memory attractor). The computation is the physical process of **energy minimization**.

---

### 14.2 The Agents: Simplification to Integrate-and-Fire

> **Summary:** For network modeling, the complex $\text{H-H}$ neuron is replaced by the simpler **Integrate-and-Fire** model, which is abstracted into a **Binary Neuron State** ($s_i \in \{+1, -1\}$). This abstraction maintains the neuron's essential **threshold decision behavior** while making large-scale network simulations computationally feasible.

#### Section Detail

The $\text{H-H}$ model, while biophysically accurate, is too intensive for large networks. The binary simplification, where $s_i$ is **active (+1)** or **silent (-1)**, is mathematically equivalent to the **Ising model**. The local update rule determines $s_i(t+1)$ based on the sign of the weighted input $\sum_j w_{ij} s_j(t)$ relative to a threshold ($\theta_i$).

#### Quiz Questions

**1. The computational abstraction of the Hodgkin–Huxley neuron used for the Hopfield Network involves replacing the continuous voltage dynamics with:**

* **A.** Continuous, stable potential $V_m$.
* **B.** **A binary state $s_i \in \{+1, -1\}$**. (**Correct**)
* **C.** Stochastic differential equations.
* **D.** The full set of gating variables.

**2. The single variable that encapsulates all the biophysical details of communication and influence between two simplified neurons, $i$ and $j$, in the Hopfield network is the:**

* **A.** Neuron state $s_i$.
* **B.** External current $I_i^{\text{ext}}$.
* **C.** **Synaptic weight $w_{ij}$**. (**Correct**)
* **D.** Membrane capacitance $C_m$.

---

#### Interview-Style Question

**Question:** The local update rule for a neuron in the Hopfield network can be seen as a sign function: $s_i(t+1) = \text{sign}(\sum_j w_{ij} s_j(t) - \theta_i)$. How does this single rule capture the concept of **Integration and Firing** from the more complex biophysical model?

**Answer Strategy:**
* **Integration:** The term $\sum_j w_{ij} s_j(t)$ represents the **weighted sum** of all incoming signals from neighboring neurons. This is the computational equivalent of the neuron **integrating** its total synaptic current.
* **Firing:** The **sign function** compares this total integrated input against the threshold ($\theta_i$). If the input exceeds the threshold, the output is positive ($\text{sign}>0 \to s_i = +1$), meaning the neuron **fires**; otherwise, it is negative ($\text{sign}\le 0 \to s_i = -1$), meaning it is silent.

---

### 14.3 The Physics Analogy: The Hopfield Network and Memory

> **Summary:** The Hopfield Network's dynamics minimize its **Energy Function**, $E(\mathbf{s}) = - \frac{1}{2} \sum_{i \neq j} w_{ij} s_i s_j + \sum_i \theta_i s_i$, which is mathematically identical to the **Ising spin glass Hamiltonian**. **Memories** ($\mathbf{s}^{(m)}$) are embedded as **low-energy attractors** by setting the symmetric weights $w_{ij}$ using the **Hebbian learning rule**.

#### Section Detail

The Hopfield Network guarantees convergence because the update rule ensures $\Delta E \le 0$ for every step. The Hebbian rule, $w_{ij} = \frac{1}{M} \sum_{m=1}^M s_i^{(m)} s_j^{(m)}$, strengthens connections between co-active neurons, literally carving out the energy valleys that store the memories. This process implements **associative recall** and **pattern completion**.

#### Quiz Questions

**1. The primary rule used to encode desired patterns $\mathbf{s}^{(m)}$ as stable memories (attractors) in the Hopfield network is known as:**

* **A.** The Boltzmann factor.
* **B.** The Integrate-and-Fire rule.
* **C.** **The Hebbian learning rule**. (**Correct**)
* **D.** The BSM equation.

**2. The dynamics of memory retrieval is guaranteed to stop changing because the energy function $E(\mathbf{s})$ is:**

* **A.** Proportional to the number of neurons, $N$.
* **B.** **Bounced from below (has a minimum value) and cannot increase with any update**. (**Correct**)
* **C.** Completely independent of the weights $w_{ij}$.
* **D.** Always zero at the stable state.

---

#### Interview-Style Question

**Question:** The memory capacity limit of the Hopfield Network is empirically defined as $M_{\text{max}} \approx 0.138 N$. Explain what happens to the energy landscape and the stored memories when a network attempts to store more than this capacity.

**Answer Strategy:** When capacity is exceeded, the **energy landscape becomes crowded**. The basins of attraction for the different stored memories begin to **overlap and interfere** with one another. This interference creates **spurious minima** (false memories) and makes the original memory attractors unstable or shallow. The network enters a **spin glass phase** where dynamics become chaotic, and recall is unreliable or results in a jumbled mixture of patterns.

---

### 14.4 The Simulation: Storing and Retrieving Patterns

> **Summary:** The simulation process involves calculating the **Hebbian weight matrix $W$** and then running the **asynchronous retrieval loop**. The asynchronous update, where a single, random neuron is updated at a time, ensures that the system performs a proper **energy descent**. The successful retrieval of a noisy input pattern demonstrates **pattern completion**.

#### Section Detail

The update rule is $s_i(t+1) = \text{sign}(\sum_{j} w_{ij} s_j(t))$, with the iteration repeating until the state stabilizes. Measuring the network **Energy $E(t)$** during retrieval confirms that the dynamics are a relaxation process, as $E(t)$ must be monotonically non-increasing. The simulation visually confirms the network's function as an **error-correcting** device.

#### Quiz Questions

**1. During the memory retrieval phase of the Hopfield Network, the dynamics are simulated using which update scheme?**

* **A.** Synchronous (all neurons update at once).
* **B.** **Asynchronous (one neuron selected and updated at a time)**. (**Correct**)
* **C.** Stochastic (Metropolis acceptance rule).
* **D.** Continuous Runge-Kutta integration.

**2. The process where the Hopfield Network is given a partial or corrupted input pattern and successfully reconstructs the full, correct stored pattern is known as:**

* **A.** Orthogonalization.
* **B.** **Pattern completion (or associative recall)**. (**Correct**)
* **C.** Critical slowing down.
* **D.** Volatility clustering.

---

#### Interview-Style Question

**Question:** The simulation requires two separate time-like processes: the initial **encoding (learning)** phase and the subsequent **retrieval (recall)** phase. Explain the key difference in the role of the network's state vector $\mathbf{s}$ between these two processes.

**Answer Strategy:**
* **Encoding (Learning Phase):** The network state $\mathbf{s}$ represents a **target memory pattern** $\mathbf{s}^{(m)}$. The network is static during this phase, and $\mathbf{s}^{(m)}$ is used to *calculate and fix* the synaptic weights $w_{ij}$ (the shape of the landscape).
* **Retrieval (Recall Phase):** The network state $\mathbf{s}$ is the **dynamical variable**. It represents the *current mental state* and is dynamically changed by the update rule. The goal is for $\mathbf{s}$ to *evolve* from a noisy starting cue until it equals a stored memory $\mathbf{s}^{(m)}$ (the final stable state).

---

## 💡 Hands-On Simulation Projects (Chapter Conclusion) 🛠️

These projects require implementing the core Hopfield Network mechanics and testing its memory properties, drawing analogies to the Ising Model.

### Project 1: Encoding and Analyzing the Weight Matrix ($W$)

* **Goal:** Implement the Hebbian learning rule and analyze the structure of the resulting weight matrix.
* **Setup:** Define three simple, orthogonal (low-overlap) binary patterns $\mathbf{s}^{(1)}, \mathbf{s}^{(2)}, \mathbf{s}^{(3)}$ on a small $N=10$ neuron network.
* **Steps:**
    1.  Implement the Hebbian learning formula $w_{ij} = \frac{1}{M} \sum_{m=1}^M s^{(m)}_i s^{(m)}_j$ to calculate the $W$ matrix.
    2.  Verify that the weight matrix is symmetric ($w_{ij}=w_{ji}$) and has zero diagonal ($w_{ii}=0$).
* ***Goal***: Establish the fundamental structure of the network's memory storage matrix.

### Project 2: Simulating Pattern Retrieval and Error Correction

* **Goal:** Demonstrate the network's ability to perform **associative recall** (memory retrieval) from a corrupted cue.
* **Setup:** Use the $W$ matrix from Project 1. Define a noisy input pattern $\mathbf{s}_{\text{input}}$ by randomly flipping $20\%$ of the bits of one stored pattern $\mathbf{s}^{(1)}$.
* **Steps:**
    1.  Implement the asynchronous retrieval loop: select a random neuron $i$, calculate $h_i = \sum w_{ij} s_j$, and update $s_i^{\text{new}} = \text{sgn}(h_i)$.
    2.  Run the simulation for 100 steps and verify that the final state matches $\mathbf{s}^{(1)}$.
* ***Goal***: Show that the final recovered state successfully corrects the $20\%$ noise in the input pattern.

### Project 3: Visualizing the Energy Landscape and Relaxation

* **Goal:** Show that the memory retrieval process is a genuine **gradient descent** (relaxation) in the energy landscape.
* **Setup:** Use the same network and retrieval simulation from Project 2.
* **Steps:**
    1.  Write a function to calculate the network's total **Energy** $E(\mathbf{s})$ at any given state $\mathbf{s}$.
    2.  During the retrieval simulation, calculate and record $E(t)$ at each time step.
* ***Goal***: Plot $E(t)$ versus time. The plot should be monotonically decreasing or constant ($\Delta E \le 0$), confirming that the asynchronous update rule always drives the system downhill toward a minimum.

### Project 4: Testing Network Capacity (The Fidelity Analogy)

* **Goal:** Demonstrate the limitation of the Hopfield Network: too many stored memories cause interference (analogous to the fidelity limits in a recording system).
* **Setup:** Use a larger $N=100$ neuron network.
* **Steps:**
    1.  Run the encoding (Hebbian rule) for a small number of patterns (e.g., $M=5$). Test retrieval fidelity (Project 2).
    2.  Run the encoding again for a large number of patterns (e.g., $M=50$). Test retrieval fidelity.
* ***Goal***: Show that recall accuracy remains high for $M=5$, but **drops significantly** for $M=50$ (since $M_{\text{max}} \approx 13$ for $N=100$), illustrating the fundamental trade-off between the number of memories and their stability.


