# 1D–1V Vlasov–Poisson Plasma Simulation

## 1. What is this project?

This project implements a **1D–1V electrostatic Vlasov–Poisson solver in Fortran** to study the evolution of a collisionless plasma using a continuous distribution function in phase space.

Unlike particle-based simulations, the Vlasov approach models the distribution of particles across position and velocity. The electric field is calculated self-consistently from the charge density using Poisson's equation.

## 2. What is it used for?

- to study the the evolution of plasma distribution functions.
- gives the understanding of self-consistent electrostatic fields and charge-density perturbations.
- exploring kinetic plasma physics and numerical simulation methods.
- analyzing electric-field energy, total energy, and particle-number conservation.
- Building a foundation for more advanced plasma simulations.

## 3. Governing equations

**Vlasov equation:**

$$
\[
\frac{\partial f}{\partial t}
+v\frac{\partial f}{\partial x}
+\frac{qE}{m}\frac{\partial f}{\partial v}=0
\]
$$ 

**Poisson's equation:**

$$
\[
\frac{\partial E}{\partial x}
=\frac{\rho}{\epsilon_0}
\]
$$

**Charge density:**

$$
\[
\rho=q(n_e-n_0), \qquad
n_e=\int f(x,v,t)\,dv
\]
$$

Here, \(f(x,v,t)\) is the distribution function, \(E\) is the electric field, \(q\) is particle charge, \(m\) is particle mass, and \(n_0\) is the uniform neutralizing-background density.

## 4. About the code

The main solver, `Vlasov_Poisson_1D1V.f90`, is written in Fortran.

- Uses a \(128\times128\) spatial–velocity grid.
- Initializes a smooth Maxwellian-like distribution.
- Evolves the distribution using explicit first-order upwind transport.
- Calculates the electrostatic field from the charge density.
- Records diagnostic quantities in `diagnostics.csv`.
- Saves the final phase-space distribution in `f_final.dat`.

Python scripts (`plot_diagnostics.py` and `plot_distribution.py`) are used to visualize the simulation outputs.

## 5. Results and analysis

The baseline simulation completed its run to \(t=10\). The reported relative mass error was approximately \(10^{-16}\), indicating excellent mass conservation for this run. The reported relative total-energy error was approximately \(5.93\times10^{-4}\), or **0.0593%**.

The electric-field energy decreased substantially during the early evolution. This may be consistent with electrostatic energy transfer and Landau-damping-like behavior, but it is not sufficient to establish Landau damping quantitatively. Further validation of the field solver, periodic boundary consistency, and grid convergence is required.

## 6. Specific assumptions and considerations

- One spatial dimension and one velocity dimension (1D–1V).
- collisionless plasma; no collision operator.
- electrostatic model; magnetic-field evolution and full Maxwell equations are not included.
- electrons evolve against a uniform, stationary positive background.
- Normalized units: \(m=1\), \(q=-1\), \(\epsilon_0=1\), and \(n_0=1\).
- Periodic spatial boundary conditions and a finite velocity domain.
- First-order upwind transport introduces numerical diffusion.
- The solver is an educational baseline and requires further numerical validation before quantitative scientific use.

## 7. Bash commands

Compile and run the Fortran solver:

```bash
gfortran Vlasov_Poisson_1D1V.f90 -o vlasov_poisson
./vlasov_poisson
```

Install Python dependencies:

```bash
python3 -m pip install numpy matplotlib
```

Generate plots:

```bash
python3 plot_diagnostics.py
python3 plot_distribution.py
```

View the output files: 

```bash
head diagnostics.csv
head f_final.dat
```

**Project status:** Educational kinetic plasma simulation under development, with Fortran-based computation and Python-based visualization.
