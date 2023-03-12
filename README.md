# Finite Element Method - Electromagnetism
Using Finite Element Method (FEM) on solving Electromagnetism problems. Some of the applications are:

- Electric Potential inside a Cylindrical Capacitor ($2b = 3.5mm$ and $2a = 1.52mm$)
- Electric Potential of a Parallel Plate Capacitor ($w = 4cm$, $t = 2mm$, $d = 1$ and $\epsilon_r = 2.2$)
- Magnetic Field ($H_z$) of TE modes in a Cylindrical Waveguide ($2a = 2cm$)
- Electric Field ($E_z$) of TM modes in a Cylindrical Waveguide ($2a = 2cm$)
- Electric Field ($E_z$) of a plane wave Scattered from Circular conductive Cylinder ($f = 300MHz$, $2a = \lambda/2$)

## Results

- Cylindrical Capacitor Potential:

![CCFEM](https://user-images.githubusercontent.com/61554467/224555540-1a3f26e4-e05d-4422-a211-3ae9c06cf1da.png)

- Parallel Plate Capacitor Potential:

![PPCFEM](https://user-images.githubusercontent.com/61554467/224555669-2e73d908-f883-447c-81bb-7255cabb625d.png)


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

For each element, we have: $$\phi = \sum_{i=1}^{3}\phi_iN_i$$
and the local integral: $$W_n = \iint_{\Omega_n}\nabla\phi'\epsilon\nabla\phi ds = \iint_{\Omega_n}(\sum_{i=1}^{3}\phi_i'\nabla N_i)\epsilon(\sum_{j=1}^{3}\phi_j\nabla N_j) ds = \sum_{i=1}^{3}\sum_{j=1}^{3}\phi_i'S_{ij}\phi_j$$

Where, $S_{ij} = \iint_{\Omega_n}\nabla N_i\epsilon\nabla N_j ds = \epsilon_n(b_ib_j + c_ic_j)A_n$ is the Local Stiffness Matric.

In our case, we don't have charge distributions in the space, so $\iint_{\Omega}\phi'\rho ds = 0$

Now, the only thing remaining is the assembly.

We seperate the nodes in groups of known and unknown potentials. For the unknowns, we create the Stifness Matrix $S_{ff}$ from the assembly of the corresponding Local Stiffness Matrices. For the known values, we create a vector $B$ that is the product of the assembly of the rest of the Local Stiffness Matrices and the known potentials.

Finally, the unknown potentials can be calculated from the linear system:
$$S_{ff}F_{f} = B$$
