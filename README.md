# Finite Element Method - Electromagnetism
Using Finite Element Method (FEM) on solving Electromagnetism problems. Some of the applications are:

- Electric Potential inside a Cylindrical Capacitor ($2b = 3.5mm$ and $2a = 1.52mm$)
- Electric Potential of a Parallel Plate Capacitor ($w = 4cm$, $t = 2mm$, $d = 1$ and $\epsilon_r = 2.2$)
- Magnetic Field ($H_z$) of TE modes in a Cylindrical Waveguide ($2a = 2cm$)
- Electric Field ($E_z$) of TM modes in a Cylindrical Waveguide ($2a = 2cm$)
- Electric Field ($E_z$) of a plane wave Scattered from Circular conductive Cylinder ($f = 300MHz$, $2a = \lambda/2$)

## Theoretical Analysis

In our problems we will use base functions to describe the quantities we want to calculate. So, their description will be in the form:
$$u(\textbf{x}) = \sum_{i = 1}^{M}u_iN_i(\textbf{x})$$
where, $u_i$ are the degrees of freedom and $N_i(\textbf{x})$ the base functions.

We will use the Simplex Coordinate System: $\zeta_i = a_i + b_ix + c_iy$ and $N_i = \zeta_i$

### Electrostatic Problem

The aim on this problem is the calculation of the electric potential in spaces with with known potentials at their boundary conditions. In our application, we will examine the case of the Cylindrical and the Parallel Plate Capacitor.

Since the nature of the problem is static, the Partial Differential Equation that we have to solve is the Poisson PDE:
$$\nabla(\epsilon\nabla\phi) + \rho = 0$$

We can restate this original problem using the weighted residual formulation (Galerkin):
$$\iint_{\Omega}\phi'(\nabla(\epsilon\nabla\phi) + \rho)ds = 0$$
where $\phi'$ are testing functions.

Using the equation: $\nabla(f\textbf{A}) = \nabla f\textbf{A} + f\nabla{\textbf{A}}$, we have:
$$\iint_{\Omega}\phi'(\nabla(\epsilon\nabla\phi) + \rho)ds = - \iint_{\Omega}\nabla\phi'\epsilon\nabla\phi ds + \iint_{\Omega} \nabla(\phi'\epsilon\nabla\phi)ds$$
and after that, with the Gauss theorem: $\iint_{\Omega}\nabla{\textbf{F}}ds = \oint_{\partial \Omega}\textbf{F}\hat{\textbf{n}}dl$ we get:
$$\iint_{\Omega}\phi'(\nabla(\epsilon\nabla\phi) + \rho)ds = - \iint_{\Omega}\nabla\phi'\epsilon\nabla\phi ds + \oint_{\partial \Omega}\phi'\epsilon\nabla\phi\hat{\textbf{n}}dl = - \iint_{\Omega}\nabla\phi'\epsilon\nabla\phi ds + \oint_{\partial \Omega}\phi'\epsilon\frac{\partial\phi}{\partial n}dl$$

In our problems due to the Neumann boundary conditions, we have $\oint_{\partial \Omega}\phi'\epsilon\frac{\partial\phi}{\partial n}dl = 0$. So, the problem takes the form:
$$\iint_{\Omega}\nabla\phi'\epsilon\nabla\phi ds = \iint_{\Omega}\phi'\rho ds$$

After that, we can switch the above equation to its discrete form (Discrete Galerkin Formulation):
$$\sum_{n=1}^{N_e}\iint_{\Omega_n}\nabla\phi'\epsilon\nabla\phi ds = \sum_{n=1}^{N_e}\iint_{\Omega_n}\phi'\rho ds$$
