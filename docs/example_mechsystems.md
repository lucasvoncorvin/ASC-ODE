# Mechanical Systems

Up to now we solved ODE systems of the form

$$
y' = f(y,t) = f(y),
$$

(for autonomous systems), coming from simple mass–spring models with one mass.

But as soon as we want to simulate many masses, springs, anchors, and constraints, the system becomes much larger.

This section explains to handle full mass–spring networks in 2D and 3D, including optional constraints.

---

# 1. Building a Full Mass–Spring Network

In a general structure (hanging chain, bridge, pendulum with a rigid bar, etc.), we have:

- multiple masses with positions $x_i(t)$  
- fixed points
- springs connecting pairs of nodes  
- optional fixed length distance constraints
- a gravity vector

All of these are stored inside:

`MassSpringSystem<dim>`

where `dim` is either `2` or `3`.

Each element (Mass, Fix, Spring, Constraint) is added using a corresponding `add(...)` method.

---

# 2. How Forces Enter the ODE

Newton’s law for mass $i$:

$$
m_i \ddot{x}_i
=
F_i^{\text{gravity}}
+ 
\sum_j F_{ij}^{\text{spring}}
+
\sum_j F_{ij}^{\text{constraint}}.
$$

## 2.1 Spring Forces

Hookean spring law:

$$
F_{ij}
=
k_{ij} \big( \|x_j - x_i\| - L_{ij} \big)
\frac{x_j - x_i}{\|x_j - x_i\|}.
$$

Gravity is constant.

## 2.2 Constraint Forces via Lagrange Multipliers

Instead of manually writing constraint forces, we enlarge the state vector:

$$
X = (x, \lambda),
$$

where $\lambda$ are constraint multipliers. 

The force vector $\vec{f}$ induced by a constraint $g(x)$ is given by:
$$
\vec{f} = \lambda \nabla g(x)
$$
where $\lambda$ is the Lagrange multiplier and $\nabla g(x)$ is the gradient of the constraint with respect to the positions.

The system of equations becomes:
$$
\begin{aligned}
F &= \nabla_x L(x, \lambda) \\
0 &= \nabla_\lambda L(x, \lambda)
\end{aligned}
$$
where $L(x, \lambda)$ is the Lagrangian including constraints.

This formulation ensures that the constraint is satisfied at all times, and the corresponding force is automatically computed.

The right-hand side is implemented in:

`MSS_Function<dim>`

which evaluates all forces, constraints, and Jacobians.

---

# 3. Time Integration: The Generalized α-Method

In this chapter we introduce the numerical method used to solve the system of equations.

Large spring networks are stiff → explicit methods blow up.  

In class, the Newmark method was introduced.  
However, for nonlinear systems it can become unstable.

The **generalized α-method** extends Newmark by adding controllable numerical damping.

### Properties:

- unconditionally stable  
- second-order accurate  
- suitable for stiff springs  
- controllable high-frequency damping  

Parameter:

- $p_\infty = 1$ → Newmark (no damping)  
- $p_\infty = 0$ → strong damping of high frequencies  

We call:

`SolveODE_Alpha(...)`

from our course ODE framework.
It expects the mass matrix (possibly constraint-modified) and the RHS function.
Both are provided by:
- `MSS_Function<dim>` — computes forces
- `ConstrainedMassMatrixFunction<dim>` — builds the mass matrix

This is wrapped inside the Python binding so the user only calls:
```py
mss.simulate(tend=10.0, steps=1000)
```
Everything else happens inside the C++ backend.

## 4.1 The Jacobian and Why We Need It

Because our system contains nonlinear springs and constraints enforced 
by Lagrange multipliers, the right-hand side

$$
F(x,\lambda)
=
\begin{pmatrix}
F_{\text{forces}}(x,\lambda) \\
C(x)
\end{pmatrix}
$$

is **nonlinear** in both $x$ and $\lambda$.  
The generalized α-method uses **implicit integration**, which means every timestep requires a 
Newton solve of the form:

$$
J \, \Delta X = -F,
$$

where

$$
J = \frac{\partial F}{\partial (x,\lambda)}.
$$

Thus we must provide the full Jacobian of the system:

$$
J =
\begin{pmatrix}
\frac{\partial F}{\partial x} & \frac{\partial F}{\partial \lambda} \\
\frac{\partial C}{\partial x} & 0
\end{pmatrix}.
$$

---

### 4.1.1 Jacobian of Springs

Each spring connecting positions $p_1$ and $p_2$ contributes a nonlinear force

$$
f = k (L - L_0) n,
\qquad
n = \frac{p_2 - p_1}{\|p_2 - p_1\|},
\qquad
L = \|p_2 - p_1\|,
$$

where $k$ is the stiffness and $L_0$ is the rest length.

The exact tangent stiffness matrix is:

$$
K = k \, n n^\top
  + \frac{k(L-L_0)}{L} \big( I - n n^\top \big),
$$

which can be decomposed into:
- **Material stiffness**: $k \, n n^\top$ (depends on stiffness)
- **Geometric stiffness**: $\frac{k(L-L_0)}{L} \big( I - n n^\top \big)$ (depends on current force)

For a spring between masses $i$ and $j$, the force on mass $i$ is $f_i = +k(L-L_0)n$ and on mass $j$ is $f_j = -k(L-L_0)n$.

The force derivatives are:

$$
\frac{\partial f_i}{\partial p_i} = -K,
\qquad
\frac{\partial f_i}{\partial p_j} = +K,
$$
$$
\frac{\partial f_j}{\partial p_i} = +K,
\qquad
\frac{\partial f_j}{\partial p_j} = -K,
$$

which enter the global Jacobian as:

- **Diagonal blocks** (e.g., $\frac{\partial f_i}{\partial p_i}$): $-K$
- **Off-diagonal blocks** (e.g., $\frac{\partial f_i}{\partial p_j}$): $+K$

Spring Jacobians give both geometric and material stiffness.

---

### 4.1.2 Jacobian of Constraints

A distance constraint enforces:

$$
C = \|p_2 - p_1\| - L_0 = 0,
$$

where $L_0$ is the fixed constraint length.

Define the unit direction vector:

$$
n = \frac{p_2 - p_1}{L},
\qquad
L = \|p_2 - p_1\|.
$$

The constraint gradient is:

$$
\frac{\partial C}{\partial p_1} = -n,
\qquad
\frac{\partial C}{\partial p_2} = +n.
$$

Constraint forces are applied via Lagrange multipliers. In the code, the gradient is stored as `dC1 = (p2-p1)/L = n`, which represents the direction from $p_1$ to $p_2$. The forces applied are:

$$
f_1 = \lambda \cdot n \quad \text{(code: } \texttt{lam * dC1}\text{)},
$$
$$
f_2 = -\lambda \cdot n \quad \text{(code: } \texttt{lam * dC2}\text{)}.
$$

The derivative of the unit direction with respect to position changes is:

$$
\frac{\partial n}{\partial (p_2 - p_1)}
= \frac{1}{L}\big(I - nn^\top\big) \equiv M.
$$

Therefore:

$$
\frac{\partial n}{\partial p_1} = -M,
\qquad
\frac{\partial n}{\partial p_2} = +M.
$$

The Jacobian block for constraint forces is:

$$
\frac{\partial f_1}{\partial p_1} = -\lambda M,
\qquad
\frac{\partial f_1}{\partial p_2} = +\lambda M,
$$
$$
\frac{\partial f_2}{\partial p_1} = +\lambda M,
\qquad
\frac{\partial f_2}{\partial p_2} = -\lambda M.
$$

In matrix form:

$$
\frac{\partial F_{\text{constraint}}}{\partial x}
=
\begin{pmatrix}
-\lambda M & +\lambda M \\
+\lambda M & -\lambda M
\end{pmatrix},
\qquad
M = \frac{1}{L}(I - n n^\top).
$$

The derivative of constraint forces w.r.t. the multiplier $\lambda$:

$$
\frac{\partial f_1}{\partial \lambda} = +n,
\qquad
\frac{\partial f_2}{\partial \lambda} = -n.
$$

In vector form:

$$
\frac{\partial F}{\partial \lambda}
=
\begin{pmatrix}
+n \\ -n
\end{pmatrix}.
$$

The constraint equation derivatives (bottom row of Jacobian):

$$
\frac{\partial C}{\partial p_1} = -n,
\qquad
\frac{\partial C}{\partial p_2} = +n,
\qquad
\frac{\partial C}{\partial \lambda} = 0.
$$

In row vector form:

$$
\frac{\partial C}{\partial x}
=
[-n^\top \;\; +n^\top],
\qquad
\frac{\partial C}{\partial \lambda} = 0.
$$

### 4.1.3 Full Constraint Jacobian Structure

For a system with $N$ masses (dimension $dN$) and $C$ constraints, the Jacobian is a $(dN+C) \times (dN+C)$ matrix.

For a single constraint $k$ between masses $i$ and $j$, combining all derivatives:

$$
J =
\begin{pmatrix}
\frac{\partial F_{\text{forces}}}{\partial x} & \frac{\partial F_{\text{forces}}}{\partial \lambda} \\
\frac{\partial C}{\partial x} & 0
\end{pmatrix}.
$$

For the specific case of **one constraint between two masses** (mass 1 and mass 2), the structure is:

$$
J =
\begin{pmatrix}
-\lambda M & +\lambda M & +n \\
+\lambda M & -\lambda M & -n \\
-n^\top     & +n^\top     & 0
\end{pmatrix},
$$


The Jacobian is implemented inside the method:

`MSS_Function<D>::evaluateDeriv(...)`

# 5.Interaction From Python
To make the system usable inside Jupyter, we use a pybind11 module (`mass_spring`).
It exposes:
- MassSpringSystem3d / MassSpringSystem2d
- Mass, Fix, Spring, Chain, DistanceConstraint
- access to lists: masses, springs, constraints 
- the simulate method

# 5.1 Python Example
Here is a minimal example:
```py

import mass_spring as ms

mss = ms.MassSpringSystem3d()
mss.gravity = [0, -9.81, 0]

# build a simple chain
f = mss.add(ms.Fix([0,0,0]))
m1 = mss.add(ms.Mass(1.0, [1,0,0]))
m2 = mss.add(ms.Mass(1.0, [2,0,0]))

mss.add(ms.Spring(1.0, 20.0, [f, m1]))
mss.add(ms.Spring(1.0, 20.0, [m1, m2]))

# optional constraint between m1 and m2
mss.add(ms.DistanceConstraint(1.0, [m1, m2]))

# run the simulation
mss.simulate(tend=5, steps=500)

```
# 5.2 Visualization in Jupyter
After simulation, the positions of the masses have been updated, and we can visualize them directly using Three.js objects in python.


Visualization in Jupyter (Short Notes)
We can use pythreejs to visualize masses as spheres and springs as line segments.
Because the C++ backend maintains up-to-date positions, animating is as simple as:

```py
for i, m in enumerate(mss.masses):
    sphere_list[i].position = m.pos
```


# 5.3 Summary of Implementation with Python

In this extension of the project we:
- generalized the mass–spring model to many nodes,
- added fixes, chains, and constraints,
- implemented the full system as a DAE with Lagrange multipliers,
- applied the α-method for stable time integration,
- wrote Python bindings so the system can be used interactively,
- created Jupyter notebooks for visualization.
- Build the module using CMake

User guideline
- Import it in Python
- Create a MassSpringSystem2d or MassSpringSystem3d
- Add masses, fixes, springs, and optionally constraints
- Call .simulate(tend, steps)
- Read out m.pos for visualization or analysis

# 6. Examples
Below is a summary of the Jupyter notebooks included with the project.

## 6.1 Double Pendulum
`mass_spring_constraint.ipynb` – Rigid Rod / Pendulum Using Constraints
This example introduces distance constraints, which enforce a fixed length between two nodes using Lagrange multipliers. We introduce this because modeling the problem with its angle becomes quite difficult for complex systems.
The constraints thereby enforces the system, s.t. The mass always stays exactly at distance l from the anchor point x_0. So it always lies on a circle with radius l, which is the rod length.


Core example
A rigid pendulum:
```py
fix = mss.add(ms.Fix([0,0,0]))
m   = mss.add(ms.Mass(1.0, [1,0,0]))

# keep them exactly 1 unit apart
mss.add(ms.DistanceConstraint(1.0, [fix, m]))

mss.simulate(5.0, 600)
```

## 6.2 Hanging Chain
`mass_spring_chain.ipynb` – Hanging Chain
This is the “intro example” notebook.
It shows the simplest nontrivial use case: a chain of masses connected by springs, hanging under gravity.
What it demonstrates:
- Creating 3D masses and fixes
- Adding a sequence of springs
- Running the α–method integrator
- Visualizing positions with pythreejs

Core idea:
We build something like:


```
Fix
 │
Spring
 │
Mass 1
 │
Spring
 │
Mass 2
 │
⋮
```


Key code pattern
```py
mss = ms.MassSpringSystem3d()
mss.gravity = [0, -9.81, 0]

fix = mss.add(ms.Fix([0,0,0]))

m_prev = fix
for i in range(10):
    m = mss.add(ms.Mass(1.0, [0, -i-1, 0]))
    mss.add(ms.Spring(1.0, 20, [m_prev, m]))
    m_prev = m

mss.simulate(5.0, 800)
```


## 6.3 Mass Spring Bridge
`mass_spring_bridge.ipynb` – A Little Cable Bridge / Lattice

Instead of a simple chain, we build a 2D grid of masses connected by springs, so the structure behaves like a flexible mini-bridge.
What it demonstrates:
- Creating structured grids of masses
- Connecting each node to neighbors
- Adding cross-bracing diagonals
- Fixing anchor nodes at both ends
- Dynamic deformation under gravity

A periodic force is applied on all nodes in z direction to simulate dynamic loads.




## 6.4 Mass Spring Kreisel
`mass_spring_kreisel.ipynb` – Spinning Spring Kreisel
This creates a kind of spinning top / rotor made of springs.
As it rotates, the spring arrangement stabilizes due to centripetal forces.

We arrange 3 masses on a circle and one fixed mass on the bottom and connect neighbors with springs.

Then we spin the structure by assigning initial velocities corresponding to a circular motion.

Key code snippet
```py
import numpy as np

N = 12
R = 1.0
masses = []

# place masses in a circle
for i in range(N):
    angle = 2*np.pi*i/N
    p = [R*np.cos(angle), R*np.sin(angle), 0]
    m = mss.add(ms.Mass(0.1, p))
    masses.append(m)

# connect neighbors
for i in range(N):
    mss.add(ms.Spring(1.0, 100.0, [masses[i], masses[(i+1)%N]]))

# give initial spinning velocity

```

# 7. Summary of Example Notebooks

| Notebook                | Description                                |
|-------------------------|--------------------------------------------|
| `mass_spring_chain`     | Basic spring networks, ropes               |
| `mass_spring_bridge`    | 2D lattices, structural stability          |
| `mass_spring_constraint`| Lagrange constraints, rigid rods           |
| `mass_spring_kreisel`   | Rotation, stiff dynamics                   |







