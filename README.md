# Quantum Simulation of Atom-Field Interactions in Dispersive Optical Lattices

This repository contains a numerical simulation of the Kapitza-Dirac effect in the adiabatic approximation, implemented in Julia. The simulation explores the quantum dynamics of atoms interacting with standing electromagnetic waves, a fundamental phenomenon in quantum optics and atomic physics.

## Physical Background

### The Kapitza-Dirac Effect

The Kapitza-Dirac effect is a quantum mechanical phenomenon where matter waves (in this case, atoms) diffract from a standing light wave. This effect is analogous to the diffraction of light by a material grating, but with the roles of light and matter reversed.

### Mathematical Model

The system is described by a set of coupled differential equations in momentum space. In the adiabatic approximation, where we consider the internal atomic dynamics to be much faster than the center-of-mass motion, the system reduces to:

$i \dot{a}_n(t) = \left( \omega_{\text{rec}} n^2 + \Omega_{\text{dim}} \right) a_n(t) + \frac{\Omega_{\text{dim}}}{2} \left( a_{n-2}(t) + a_{n+2}(t) \right)$

where:

- $a_n(t)$  represents the probability amplitude for the atom to have momentum $n\hbar k$
- $ω_{rec}$  is the recoil frequency
- $Ω_{dim}$ is the effective Rabi frequency, given by  $Ω_{dim} = Ω²/(2δ)$

## Implementation

### Dependencies

- `DifferentialEquations.jl`: For solving the system of ODEs
- `ParameterizedFunctions.jl`: For defining the ODEs symbolically
- `Plots.jl`: For visualization

### Key Components

1. **System Definition**
   
   ```julia
   my_ode10 = @ode_def begin
       # System of 21 coupled differential equations
       # representing momentum states from -10 to +10
       ...
   end Ω ωᵣ
   ```

2. **Solution Extraction**
   
   ```julia
   function extract_solutions(solutions)
       # Efficiently extracts probability amplitudes
       # from the ODE solution
       ...
   end
   ```

3. **Probability Calculation**
   
   ```julia
   function calculate_probabilities(solutions)
       # Calculates transition probabilities
       # |a_n(t)|² for each momentum state
       ...
   end
   ```

## Usage

1. Set the physical parameters:
   
   ```julia
   const ω = 2π * 73.7e3  # Frequency for Li-6
   const V₀ = 50ω        # Potential strength
   ```

2. Define initial conditions:
   
   ```julia
   u₀ = Complex{Float64}[1; zeros(20)]  # Initial state in ground state
   ```

3. Run the simulation:
   
   ```julia
   prob = ODEProblem(my_ode10, u₀, tspan, p)
   sol = solve(prob, abstol=1e-15, reltol=1e-20)
   ```

4. Visualize results:
   
   ```julia
   solutions = extract_solutions(sol.u)
   probs = calculate_probabilities(solutions)
   plot_probabilities(sol.t, probs)
   ```

## Results

The simulation produces plots showing the time evolution of transition probabilities between different momentum states. These demonstrate the characteristic oscillatory behavior of atoms diffracting from the standing wave.

## Physical Parameters

- For Lithium-6 (⁶Li):
  - Recoil frequency: $ωᵣ = 2π * 73.7$ kHz
  - Potential strength: $V₀ = 50ωᵣ$
- For Rubidium-87 (⁸⁷Rb):
  - Recoil frequency: $ωᵣ = 2π * 3.8$ kHz

## References

1. Báez de la Luz, N. (2023). Interacción átomo-campo en redes ópticas dispersivas [Bachelor's thesis, Universidad Nacional Autónoma de México]. UNAM Repository. [TESIUNAM - Vista completa del registro](https://acortar.link/QdmYmz)
2. Kapitza, P. L., & Dirac, P. A. M. (1933). The reflection of electrons from standing light waves. Mathematical Proceedings of the Cambridge Philosophical Society, 29(2), 297-300.

## Future Work

Potential extensions of this project include:

- Implementation of non-adiabatic effects
- Addition of decoherence and dissipation
- Extension to multiple atomic levels
- Integration with machine learning for parameter optimization
- Development of a web-based interactive simulator

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
