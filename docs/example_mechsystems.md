# Mass–Spring Systems as ODEs and DAEs

Up to now we solved ODE systems of the form

$$
y' = f(y,t) = f(y),
$$

(for autonomous systems), coming from simple mass–spring models with one mass.

But as soon as we want to simulate many masses, springs, anchors, and constraints, the system becomes much larger.

This section explains **how we extended the toolbox** to handle full mass–spring networks in 2D and 3D, including optional constraints.

---

# 1. Building a Full Mass–Spring Network

In a general structure (hanging chain, bridge, pendulum with a rigid bar, etc.), we have:

- multiple masses with positions $x_i(t)$  
- fixes (immovable anchor points)  
- springs connecting pairs of nodes  
- optional fixed-length **distance constraints**  
- a global gravity vector

All of these are stored inside:

MassSpringSystem<dim>

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

Spring forces follow Hooke’s law:

$$
F_{ij}
=
k_{ij} \big( \|x_j - x_i\| - L_{ij} \big)
\frac{x_j - x_i}{\|x_j - x_i\|}.
$$

Gravity is constant.

What is new is that **constraints** are handled using **Lagrange multipliers**.

We enlarge the state vector:

$$
X = (x, \lambda),
$$

where $\lambda$ are constraint multipliers.

This turns the dynamics into a **Differential–Algebraic Equation (DAE)**.

The right-hand side is implemented in:

MSS_Function<dim>

which evaluates all forces, constraints, and Jacobians.

---

# 3. The Extended State Vector

For $N$ masses and $C$ constraints the sizes are:

- positions: $dN$  
- velocities: $dN$  
- accelerations: $dN$  
- multipliers: $C$

Thus the state stored in the solver is:
x[0 : dN] → mass positions
x[dN : dN + C] → constraint multipliers λ

`getState` and `setState` convert this automatically.

---

# 4. Solving the System: The Generalized α-Method

For large, stiff networks (which is exactly what happens with many springs), explicit methods blow up quickly. In class the Newmark method was introduced. Since it can lead to instabilities for non linear functions, an extension of this method leads to the generalized alpha method to allow controlled numerical damping, where the parameter pinfinity specifies the damping for high frequency behavior. If we set it to 1, we recover the Newmark method ( no numerical damping) and if we set it to 0, we achieve maximum high frequency damping.
So in this case, we switched to the generalized α-method, an implicit integrator that is:

- unconditional stability  
- tunable high-frequency damping  
- second-order accuracy  
- robustness for stiff springs  

Parameter $p_\infty$ controls the damping:

We use the solver function:

SolveODE_Alpha(...)

from our course ODE framework.
It expects the mass matrix (possibly constraint-modified) and the RHS function.
Both are provided by:
- MSS_Function<dim> — computes forces
- ConstrainedMassMatrixFunction<dim> — builds the mass matrix

This is wrapped inside the Python binding so the user only calls:

mss.simulate(tend=10.0, steps=1000)

Everything else happens inside the C++ backend.
Interaction From Python
To make the system usable inside Jupyter, we use a pybind11 module (mass_spring).
It exposes:
- MassSpringSystem3d / MassSpringSystem2d
- Mass, Fix, Spring, Chain, DistanceConstraint
- access to lists: masses, springs, constraints 
- the simulate method


## 5.3 The Jacobian and Why We Need It

Because our system contains nonlinear springs and holonomic constraints enforced 
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

Thus we must provide the full Jacobian of the DAE system:

$$
J =
\begin{pmatrix}
\frac{\partial F}{\partial x} & \frac{\partial F}{\partial \lambda} \\
\frac{\partial C}{\partial x} & 0
\end{pmatrix}.
$$

---

### 5.3.1 Jacobian of Springs

Each spring contributes a nonlinear force

$$
f = k (L - L_0) n,
\qquad
n = \frac{p_2 - p_1}{\|p_2 - p_1\|}.
$$

The exact tangent stiffness matrix is:

$$
K = k n n^\top
  + \frac{k(L-L_0)}{L} \big( I - n n^\top \big),
$$

which enters the Jacobian block:

- $+\;K$ on the diagonal of each connected mass  
- $-\;K$ on the off-diagonal blocks  

Spring Jacobians give geometric + material stiffness.

---

### 5.3.2 Jacobian of Constraints

A distance constraint has:

$$
C = \|p_2 - p_1\| - L_0,
\qquad
n = \frac{p_2 - p_1}{L}.
$$

Its gradient is:

$$
\frac{\partial C}{\partial p_1} = n,
\qquad
\frac{\partial C}{\partial p_2} = -n.
$$

Constraint forces are:

$$
f_1 = \lambda n,
\qquad
f_2 = -\lambda n.
$$

The derivative of the unit direction is:

$$
\frac{\partial n}{\partial p}
= \frac{1}{L}\big(I - nn^\top\big).
$$

So the Jacobian block for constraint forces is:

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

Derivative w.r.t. $\lambda$:

$$
\frac{\partial F}{\partial \lambda}
=
\begin{pmatrix}
n \\ -n
\end{pmatrix}.
$$

Bottom row:

$$
\frac{\partial C}{\partial x}
=
[n \;\; -n],
\qquad
\frac{\partial C}{\partial \lambda} = 0.
$$

## 6. Full Constraint Jacobian

Combining:

- the force derivatives w.r.t. $x$
- force derivatives w.r.t. $\lambda$
- constraint derivatives w.r.t. $x$
- and the zero block for $\partial C / \partial \lambda$

we obtain the exact Newton Jacobian:

$$
J =
\begin{pmatrix}
-\lambda M & +\lambda M & n \\
+\lambda M & -\lambda M & -n \\
n^\top     & -n^\top     & 0
\end{pmatrix}.
$$

This is the complete analytic Jacobian for one distance constraint.

In a multi-mass system, these blocks are placed into the appropriate rows/columns corresponding to the constrained masses.

---

### 5.3.3 Where the Jacobian Is Implemented in the Code

The Jacobian is implemented inside the method:

MSS_Function<D>::evaluateDeriv(...)

Here is a minimal example:

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

After simulation, the positions of the masses have been updated, and we can visualize them directly using Three.js objects in python.


Visualization in Jupyter (Short Notes)
We can use pythreejs to visualize masses as spheres and springs as line segments.
Because the C++ backend maintains up-to-date positions, animating is as simple as:

for i, m in enumerate(mss.masses):
    sphere_list[i].position = m.pos


This allowed us to build hanging chains, constrained pendulums, a little “bridge”, and a spinning spring–kreisel model.

Summary of This Extension
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

Examples

mass_spring_chain.ipynb – Hanging Chain
This is the “intro example” notebook.
It shows the simplest nontrivial use case: a chain of masses connected by springs, hanging under gravity.
What it demonstrates:
- Creating 3D masses and fixes
- Adding a sequence of springs
- Running the α–method integrator
- Visualizing positions with pythreejs

Core idea:
We build something like:

Fix
 |
Spring
 |
Mass1
 |
Spring
|
Mass2
 ...


Key code pattern
mss = ms.MassSpringSystem3d()
mss.gravity = [0, -9.81, 0]

fix = mss.add(ms.Fix([0,0,0]))

m_prev = fix
for i in range(10):
    m = mss.add(ms.Mass(1.0, [0, -i-1, 0]))
    mss.add(ms.Spring(1.0, 20, [m_prev, m]))
    m_prev = m

mss.simulate(5.0, 800)


You see the chain sag and oscillate — a sort of sanity check that the system behaves as a physical rope.


10.2 mass_spring_bridge.ipynb – A Little Cable Bridge / Lattice

Instead of a simple chain, we build a 2D grid of masses connected by springs, so the structure behaves like a flexible mini-bridge.
What it demonstrates:
- Creating structured grids of masses
- Connecting each node to neighbors
- Adding cross-bracing diagonals
- Fixing anchor nodes at both ends
- Dynamic deformation under gravity

Structure

We create a W×H lattice:
●───●───●───●   ← fixed on left + right
│  /│  /│  /
●───●───●───●


Key code pattern
W, H = 6, 2
nodes = []

for j in range(H):
    row = []
    for i in range(W):
        m = mss.add(ms.Mass(0.2, [i, -j, 0]))
        row.append(m)
    nodes.append(row)

# fix the left and right ends
mss.add(ms.Fix([0,0,0]))
mss.add(ms.Fix([W-1,0,0]))

# connect horizontal, vertical, and diagonal springs

When simulated, the bridge bends and vibrates realistically — and this already shows why stable implicit integration is essential.

10.3 mass_spring_constraint.ipynb – Rigid Rod / Pendulum Using Constraints
This example introduces distance constraints, which enforce a fixed length between two nodes using Lagrange multipliers. We introduce this because modeling the problem with its angle becomes quite difficult for complex systems.
The constraints thereby enforces the system, s.t. The mass always stays exactly at distance l from the anchor point x_0. So it always lies on a circle with radius l, which is the rod length.

What it demonstrates:
- Adding DistanceConstraint objects
- Mixing springs and constraints
- Understanding how constraints give rigid-body-like behavior


Core example
A rigid pendulum:
fix = mss.add(ms.Fix([0,0,0]))
m   = mss.add(ms.Mass(1.0, [1,0,0]))

# keep them exactly 1 unit apart
mss.add(ms.DistanceConstraint(1.0, [fix, m]))

mss.simulate(5.0, 600)

The result is a clean swing, not a stretchy spring.
This notebook demonstrates how constraints do not require any user-side derivative computations — the library handles the DAE system internally.


10.4 mass_spring_kreisel.ipynb – Spinning Spring Kreisel
This creates a kind of spinning top / rotor made of springs.
As it rotates, the spring arrangement stabilizes due to centripetal forces.

What it demonstrates:
- Building circular spring structures
- Understanding rotation stabilization
- The effect of stiffness vs. rotation speed
- How our solver handles stiff, fast-moving systems

Geometry
We arrange masses on a circle and connect neighbors with springs:
    ●───●
   /       \
  ●         ●
   \       /
     ●───●

Then spin the structure by assigning initial velocities.

Key code snippet
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

The kreisel stabilizes, showing off the robustness of the α-method.


10.5 double_chain_pendulum.ipynb – Two Pendulums Linked in a Chain
This notebook combines everything:
- multiple masses
- constraints
- springs
- gravity
- chaotic motion

The system is a “double pendulum” but built using our API rather than closed-form equations.
This makes the motion unpredictable and visually interesting.

What it demonstrates:
- Chaining constraints to make articulated structures
- Nonlinear, chaotic dynamics
- Importance of small time steps

Core idea
Fix  ●
      \
       (constraint)
        \
         ● Mass 1
          \
           (constraint)
            \
             ● Mass 2

Key snippet
f = mss.add(ms.Fix([0,0,0]))
m1 = mss.add(ms.Mass(1.0, [1,0,0]))
m2 = mss.add(ms.Mass(1.0, [2,0,0]))

mss.add(ms.DistanceConstraint(1.0, [f, m1]))
mss.add(ms.DistanceConstraint(1.0, [m1, m2]))

mss.simulate(10.0, 2000)

The motion looks chaotic — just like a theoretical double pendulum should — but here it comes purely from the mass–spring–constraint model.

11. Summary: What Users Can Learn from Our Examples
Each notebook highlights a different aspect of the library:
Notebook                           Focus
mass_spring_chain                  Basic spring networks, ropes

mass_spring_bridge                 2D lattices, structural stability

mass_spring_constraint             Lagrange constraints, rigid rods

mass_spring_kreisel                Rotation, stiff dynamics

double_chain_pendulum              Chains of constraints, chaotic motion






