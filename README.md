
# Simulating Complex Systems

[![Documentation](https://img.shields.io/badge/docs-live-brightgreen)](https://bigbookofcomputing.github.io)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
[![MkDocs](https://img.shields.io/badge/Built%20with-MkDocs-blue)](https://www.mkdocs.org/)

> **Volume II** of the *Big Book of Computing* series

## 📖 About

**Simulating Complex Systems** is a comprehensive guide to modeling and simulating complex phenomena across physics, finance, and biology. This volume addresses the "many-body problem"—the central challenge in understanding systems where interactions between components give rise to emergent collective behavior that cannot be predicted from individual elements alone.

From Monte Carlo methods to agent-based models, from molecular dynamics to stochastic calculus, this book provides a unified computational toolkit for exploring complex systems where traditional analytical approaches fail.

## 🌟 Why This Book?

Real-world systems—whether spin lattices, financial markets, or cellular networks—are too complex for closed-form solutions. Mean-field approximations break down when fluctuations, correlations, and nonlinear interactions dominate. This book teaches you how to:

- **Simulate** systems computationally when analytical solutions don't exist
- **Sample** high-dimensional spaces efficiently using Monte Carlo methods
- **Model** emergent behavior from simple local rules
- **Analyze** stochastic and deterministic dynamics across scales
- **Bridge** physics, finance, and biology through unified computational frameworks

## 🎯 What's Inside

This book is organized into three major parts, each focusing on a different computational paradigm:

### Part I: Stochastic Simulation & Statistical Mechanics

Monte Carlo methods and importance sampling for systems in equilibrium.

- **Chapter 1**: Foundations of Stochastic Simulation — From simple sampling to Metropolis-Hastings
- **Chapter 2**: The Ising Model — Phase transitions and spontaneous symmetry breaking
- **Chapter 3**: Lattice Gauge Theory — Simulating quantum field theories non-perturbatively
- **Chapter 4**: Monte Carlo Option Pricing — Pricing exotic derivatives and path-dependent options
- **Chapter 5**: Stochastic Systems Biology — The Gillespie algorithm and transcriptional bursting
- **Chapter 6**: Advanced Monte Carlo Methods — Cluster algorithms, parallel tempering, and Wang-Landau

### Part II: Deterministic & Stochastic Dynamics

Time evolution of systems through differential equations and stochastic processes.

- **Chapter 7**: Molecular Dynamics — Classical motion of interacting atoms with Verlet integration
- **Chapter 8**: The Stochastic Calculus — SDEs, Wiener processes, and Itō's lemma
- **Chapter 9**: Black-Scholes-Merton — From stochastic paths to deterministic pricing PDEs
- **Chapter 10**: Computational Neuroscience — Hodgkin-Huxley equations and action potentials

### Part III: Agent-Based & Network Models

Bottom-up modeling of emergent behavior through interacting agents.

- **Chapter 11**: Agent-Based Model Framework — Emergence from simple local rules
- **Chapter 12**: Agent-Based Market Models — Herding, bubbles, and fat-tailed distributions
- **Chapter 13**: Collective Behavior & Pattern Formation — Turing patterns and reaction-diffusion
- **Chapter 14**: Computational Neuroscience II — Neural networks, attractors, and associative memory
- **Chapter 15**: Graph and Network Models — Topology, hubs, and network robustness

## 🚀 Getting Started

### View the Book Online

The complete book is available online at: **[https://bigbookofcomputing.github.io](https://bigbookofcomputing.github.io)**

### Build Locally

To build and serve the documentation locally:

1. **Clone the repository**
   ```bash
   git clone https://github.com/bigbookofcomputing/simulation.git
   cd simulation
   ```

2. **Install dependencies**
   ```bash
   pip install mkdocs-material
   pip install mkdocs-minify-plugin
   ```

3. **Serve locally**
   ```bash
   mkdocs serve
   ```
   
   Then open your browser to `http://127.0.0.1:8000`

4. **Build static site**
   ```bash
   mkdocs build
   ```

5. **Deploy to GitHub Pages**
   ```bash
   mkdocs gh-deploy
   ```

## 📚 Book Structure

Each chapter follows a three-part learning framework:

- **📖 Essay** — Deep conceptual understanding with theory and context
- **📘 WorkBook** — Hands-on exercises to build intuition through practice
- **💻 CodeBook** — Runnable implementations with detailed comments

This structure ensures you don't just read about simulation—you actively build, test, and explore complex systems.

## 🔬 Key Themes

This book emphasizes several unifying concepts across domains:

### The Physics-Finance-Biology Connection
- **Ising Model** ↔️ Market herding ↔️ Cell signaling
- **Partition functions** ↔️ Option pricing ↔️ Protein folding
- **Phase transitions** ↔️ Market crashes ↔️ Pattern formation

### Computational Paradigms
- **Equilibrium sampling** — When the goal is to characterize distributions (Monte Carlo)
- **Dynamical evolution** — When the goal is to simulate time evolution (MD, SDEs)
- **Emergent behavior** — When the goal is to discover collective phenomena (ABMs)

### Practical Skills
- Dealing with **thermalization** and **autocorrelation**
- Implementing **periodic boundary conditions**
- Analyzing **convergence** and **ergodicity**
- Measuring **observables** from noisy data
- Optimizing with **neighbor lists** and **variance reduction**

## 🛠️ Technologies Used

- **MkDocs Material** — Modern documentation framework
- **MathJax** — Mathematical typesetting
- **Python** — Primary language for implementations
- **NumPy/SciPy** — Numerical computing
- **Matplotlib** — Visualization

## 🎓 Who Should Read This Book?

This book is designed for:

- **Physics students** learning computational and statistical mechanics
- **Quantitative finance practitioners** modeling markets and pricing derivatives
- **Computational biologists** simulating cellular and neural systems
- **Data scientists** interested in stochastic processes and complex systems
- **Anyone** curious about how computation reveals emergent phenomena

**Prerequisites**: Familiarity with calculus, basic probability, and programming. Volume I (Foundation of Computational Science) provides helpful background on numerical methods.

## 🤝 Contributing

We welcome contributions! Whether it's:

- Fixing errors or improving explanations
- Adding new examples or case studies
- Suggesting additional topics
- Reporting issues

Please feel free to open an issue or submit a pull request.

## 📄 License

This work is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## 🌟 About the Big Book of Computing

This is **Volume II** of the *Big Book of Computing* series:

- **Volume I**: [Foundation of Computational Science](https://github.com/bigbookofcomputing/foundation) — Numerical methods and computational foundations
- **Volume II**: **Simulating Complex Systems** — Monte Carlo, dynamics, and agent-based models (this volume)
- **Volume III**: Optimization & Machine Learning *(coming soon)*
- **Volume IV**: Quantum Computing *(coming soon)*

## 📧 Contact

- **Website**: [https://bigbookofcomputing.github.io](https://bigbookofcomputing.github.io)
- **GitHub**: [https://github.com/bigbookofcomputing](https://github.com/bigbookofcomputing)
- **Twitter**: [@bigbookofcomputing](https://x.com/bigbookofcomputing)

---

**Built with ❤️ for understanding emergence, randomness, and complexity**
