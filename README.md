# Dual Adaptation Random Matrix Theory

This is a theoretical neuroscience project investigating two forms of adaptation found in the brain: (1) weight adaptation and (2) bias adaptation. Mechanistically these most closely align with short-term synaptic depression and spike frequency adaptation, respectively. We model the effect of dual adaptation on continuous time recurrent neural networks.

## Project Overview

Our neural network model includes:
- **Excitatory and inhibitory neurons** obeying Dale's principle
- **Rectifying (ReLU) neurons** with nonlinear dynamics
- **Adaptation only in excitatory neurons** - when there is high activity, weight and bias adaptation reduce the output of the E neurons, allowing the I neurons to "win" (i.e., dynamically balance the E neurons)

## Research Questions

Our primary research focus examines:

Motivating question: How are global dynamical properties maintained through **biologically plausible** (i.e., local) development rules? The stability of LTI dynamical systems (i.e, systems of first-order differential equations) are 100% determined by connectivity (i.e., eigenvalues). Is this true for biomimetic recurrent nonlinear networks?

Hypothesis: We hypothesize that dual adaptation helps decouple global dynamical properties from the details of connectivity. We expect that the time constants of adaptation will play a role.

## Methodology

We investigate the worst case scenario where synaptic weights are drawn from a random distribution, an approach used in **random matrix theory**. We extend concepts of random matrix theory to our nonlinear RNN by:
- Simulating random networks
- Computing the Lyapunov spectrum
- Determining the level of chaos with the **Kolmogorov-Sinai entropy**

## Current Status and Roadmap

### ✅ Completed
- [x] Created a **Lorenz test system** with known largest Lyapunov exponent for algorithm verification
- [x] Implemented **Benettin's method** for Lyapunov exponent calculation

### 🔄 In Progress / Planned
- [ ] Implement **QR factorization** to compute the full Lyapunov spectrum
- [ ] Update README with **LaTeX equations** for the RNN model
- [ ] Write **MATLAB function** for the nonlinear RNN
- [ ] Create **example implementations** using the RNN function
- [ ] Perform **sensitivity analysis** of connectivity distribution parameters
- [ ] Incorporate figures into manuscript (already written based on prior code base)

## License

MIT License 