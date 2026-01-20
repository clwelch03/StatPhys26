# Brownian Motion as a Free Energy Minimizing Process

Consider a suspension of a large number of ideal beads:

![Figure](https://cdn.mathpix.com/cropped/2025_03_13_36c7a3b802762c059e2bg-01.jpg?height=774&width=1398&top_left_y=802&top_left_x=394)

$$
\xi = 6 \pi \eta_{3} \sigma
$$

where $\xi$ is the friction coefficient.

The particle density at $x$ and time $t$ is denoted as $c(x, t)$.

### Goal:
Construct a phenomenological model for the dynamics:

$$
C(x, t) \longrightarrow C_{eq}(t) \sim e^{-\frac{U(x)}{k_B T}}.
$$

We will see that Brownian dynamics continuously lower the Free Energy.



*Conservation of particle number* demands that the dynamical equation of motion is described by the continuity equation

$$
\frac{\partial c(x,t)}{\partial t} = -\frac{\partial j(x,t)}{\partial x} \;
$$ (continuity-eq)

The probability current $j(x, t)$ (rate of particles passing through $x$ per unit time) results from:

1. **Deterministic motion:**
   $$
   j_{\text{det}} = c \cdot \bar{v} = -\frac{c}{\xi} \frac{\partial U}{\partial x}
   $$

2. **Random motion:**
   $$
   j_{\text{diff}} = -D \frac{\partial c}{\partial x}
   $$
   (Fick’s Law)

![Figure](https://cdn.mathpix.com/cropped/2025_03_13_36c7a3b802762c059e2bg-02.jpg?height=378&width=1253&top_left_y=1477&top_left_x=589)

Both contributions can be summarized by the equation

$$
j = c \cdot v_{fl} = 0
$$

in terms of a net velocity

$$
v_{fl} = -\xi^{-1} \frac{\partial}{\partial x} \left( U + D \xi \ln c \right) = -\xi^{-1} \frac{\partial \mu}{\partial x}
$$

where $\mu(x)$ is the local chemical potential of non-interacting particles:

$$
\mu(x) = k_B T \ln c + U(x).
$$

Inserting our expression for the current into {eq}`continuity-eq`  leads to the generalized diffusion equation:

$$
\frac{\partial c}{\partial t} =\frac{\partial}{\partial x}\left(c \frac{\partial_x \mu}{\xi}\right)= \frac{\partial}{\partial x} \left( \frac{1}{\xi} \left( D \xi \frac{\partial c}{\partial x} + C \frac{\partial U}{\partial x} \right) \right).
$$


### Equilibrium $t \to \infty$:


$$
\frac{\partial c}{\partial t} = -\partial_x j \to 0
$$

which implies that the flux $j$ is constant and must be zero. Thus,

$$
\mu = U + D \xi \ln c = \text{const}.
$$

From this, the equilibrium distribution follows as

$$
C_{eq} \propto \exp \left( -\frac{U(x)}{D \xi} \right).
$$

Since we have to demand that this is identical to the the Boltzmann distribution, 

$$
C_{eq} \propto \exp \left( -\frac{U(x)}{k_B T} \right)\;,
$$

we have to demand

$$
D=\frac{k_B T}{\xi} \;.
$$

This relationship is called the Stokes-Einstein relation. It is another example of a fluctuation-dissipation relationship.

To aid the interpretation of $D$, we'll anticipate that the mean squared displacement of a particle diffusing for a time $t$ is given by:

$$
\langle \Delta x^2 (t) \rangle = 2 D t.
$$

Examples: For a spherical particle of radius $R$, the friction coefficient is:

$$
\xi = 6 \pi \eta R.
$$

For typical proteins, $R \approx 2$ nm and $D \approx 100 \mu m^2/s$, while for bacteria, $D_{bac} \approx 0.5 \mu m^2/s$.

### Note

We found that the current $j$ is proportional to the gradient of our good old chemical potential, $j\propto \nabla \mu$. This automatically implies $j=0$ when $\mu=$const., which is our condition for thermodynamic equilibrium. 

This suggests a generalization of the generalized diffusion equation. Introduce a probability distribution function $\psi(\{x_i\}_i)$ of a set of degrees of freedom $\{x_i\}_{i=1\dots N}$, which could be, e.g., the monomer positions of a polymer or membrane, etc..

Since the probability is locally conserved, we demand a continuity equation

$$
\partial_t\Psi(\{x_i\}_i) =-\sum_i\partial_i j_i \;,
$$

where $j_i$ is the probability current in the $i$ direction. It may be expressed as

$$
j_i = \Psi v_i
$$

in terms of an effective velocity $v_i$,

$$
v_i = \sum_j H_{ij} F_j
$$

where $H_{ij}$ is the so-called mobility and $F_j$ is a force acting on the $j^{th}$ d.o.f.,

$$
F_j=-\partial_j \mu
$$

which depends on gradients of the chemical potential
 
$$
\mu=U + k_B T \ln(\Psi)
$$






This equation is also called Smolukowski equation.



![Figure](https://cdn.mathpix.com/cropped/2025_03_13_36c7a3b802762c059e2bg-05.jpg?height=325&width=785&top_left_y=1513&top_left_x=605)


## Free Energy Minimization

Define the functional:

$$
F[c] = \int dx \, c(x, t) \mu(x).
$$

The free energy is monotonically decreasing:

$$
\frac{dF}{dt} = \int dx \left[ \mu \frac{\partial c}{\partial t} + c \frac{\partial \mu}{\partial t} \right].
$$

Using integration by parts and boundary conditions at infinity, we obtain:

$$
\frac{dF}{dt} = -\int dx \left[ c (\partial_x\mu)^2 + k_B T \partial_t c \right] \;,
$$

which is negative unless $\mu=$const., i.e., we have reached equilibrium.

N.B.: The dynamics is time irreversible because we have coarse-grained dynamics that is chaotic at the microscopic scales. (see next lecture on Ergodicity).

### Basic solution of the diffusion equation:

$$
\begin{equation*}
U=0 \Rightarrow \frac{\partial c\left(x, t\right)}{\partial t}=D \partial_{x}^{2} c\left(x, t \right) \tag{1}
\end{equation*}
$$ (diff-equ)

![](https://cdn.mathpix.com/cropped/2025_03_18_db9d17bf3c7ddba512abg-5.jpg?height=896&width=1304&top_left_y=995&top_left_x=307)

Assume point-like initial conditions, $c(x,0)=\delta(x)$.

Fourier transform of equation {eq}`diff-equ`:
$$
\begin{aligned}
& \quad c(x, t)=\int_{-\infty}^{\infty} \frac{d q}{2 \pi} e^{-i q x} c(q, t) \\
& \frac{\partial c(q, t)}{\partial t}=-D q^{2} c(q, t) .
\end{aligned}
$$ (diff-equ-ft)

Now, this is an ordinary differential equation (in contrast to the pdf {eq}`diff-equ`).
The solution is
$$
c(q, t)=c(q, 0) e^{-D_{q}^{2} t}
$$


$c(q,0)$ follows from Fourier transforming the delta function initial condition:
$$
\begin{aligned}
& c(q, 0)=\int_{-\infty}^{\infty} d x e^{i q x} c(x, t) \\
&=\int_{-\infty}^{\infty} d x e^{i q x} \delta(x)=1 \\
& \Rightarrow c(q, t)=e^{-D q^{2} t}
\end{aligned}
$$

Fourier back transform:
$$
\begin{aligned}
c(x, t) & =\int_{-\infty}^{\infty} \frac{d q}{2 \pi} c(q, t) e^{-i q x}= \\
& =\int_{-\infty}^{\infty} \frac{d q}{2 \pi} e^{-D q^{2} t-i q x} \\
& =\int \frac{d q}{2 \pi} e^{-(\underbrace{\sqrt{D T} q-\frac{i x}{2 \sqrt{D t}}}_{\equiv Q})^{2}-\frac{x^{2}}{4 D t}} \\
& =e^{-\frac{x^{2}}{4 D t}} \frac{1}{2 \pi \sqrt{D t}} \underbrace{\int_{-\infty}^{\infty} d Q e^{-Q^{2}}}_{\sqrt{\pi}}
\end{aligned}
$$

$$
\Rightarrow \begin{align*}
& c(x, t)=\frac{e^{-\frac{x^{2}}{4 D t}}}{\sqrt{4 \pi D t}}=\frac{e^{-\frac{x^{2}}{2 \sigma_{t}^{2}}}}{\sqrt{2 \pi \sigma_{t}^{2}}}  \\
& \text { with } \sigma_{t}^{2}=\left\langle(x(t)-x(\theta))^{2}\right\rangle=2 D t
\end{align*}
$$ (greens-fct)

N. B.: This is a Green's function solution. The solution to any initial condition can be fond by a convolution with the Greensfunction. 

N.B:: In our examples (sofar), all particles are non-interacting.
$$
\Longrightarrow c(x, t)=N \cdot p(x, t)
$$
with $p(x, t)$ being the probability of a single particle to be at $x, t$.

#### Let's put in numbers.
i) Nerve Cell

![](https://cdn.mathpix.com/cropped/2025_03_18_5e6b0bbe0ca706fbe34cg-1.jpg?height=356&width=1053&top_left_y=797&top_left_x=526)

Is diffusion fast enough to transport vesicles?

$$
\begin{aligned}
D & =0,5 \mu^{2} / s=\frac{k T}{6 \pi y R_{k}} ; \\
\frac{\left\langle x^{2}\right\rangle}{2 D} & =t=12 \text { days! }
\end{aligned}
$$
$\Rightarrow$ Nerve cells need active transport.

ii) $E.coli$ has a length of about $L\approx 1\mu$m. So, the time to diffuse across the length of the bacterium is $t\approx L^2 /(2D)$. Large molecules inside the cells have nm radii, e.g., the RNA-polymerase has a radius $R\approx 2.5$ nm, and hence a diffusivity of $D\approx 100 \mu m^2$/s. The diffusion time thus is $t\approx L^2 /(2D)\approx 2$ms. This time is much shorter than the 20min life cycle of $E. coli$. So, diffusion is an effective transport mechanism for bacteria.


### Applications:

#### Diffusion sets fundamental limits to Rates of biological aud chemical reactions.

![](https://cdn.mathpix.com/cropped/2025_03_18_5e6b0bbe0ca706fbe34cg-2.jpg?height=865&width=1297&top_left_y=676&top_left_x=625)

How many ligands are captured per sec?

Use diffusion equation to model the ligands:
$$
\partial_{t} c=D \nabla^{2} c=D\left(\partial_{x}^{2} c+\partial_{y}^{2} c+\partial_{z}^{2} c\right)
$$

$$
\partial_{t} c=D \frac{\psi}{r^{2}} \partial_{r}\left(r^{2} \partial_{r} c\right) ; c=c\left(r, t\right)
$$

Steady state: $0=\partial_{t} c=0 \quad \partial_{r} c_{s}=\frac{\text { cost }}{r^{2}}$ for $r>R$

$$
c_{s}(r)=\frac{-A}{r}+B=-\frac{A}{r}+c_{0} \Rightarrow c_{s}(\infty)=c_0
$$

i) Perfect absorber limit: $c_{s}(R)=0$.

$$
\Rightarrow c_{s}(r)=c_{0}\left(1-\frac{R}{r}\right) .
$$

particle current $j(R)=|\vec{j}(\vec{R})|=D \partial_{r} c(r, t)=\frac{D c_0 R}{R^{2}}$.

1D: $j(x,t)=-D\partial_x c$

3D: $j(x,t)=-D \nabla c(\vec r,t)$

Total flux of captured particles is

$$
\Phi = j(R) 4\pi R^2=4\pi D c_0 R
$$

(Check dimensions!)

#### Universal bound on chemical reaction rates

Key question: How fast is $A+B \xrightarrow{k} C$ if it is limited by $A-B$ collision rates?


In the extreme case, where each A-B collision leads to a C particle we have $k=c_{A}\tilde{k}_{A-B}$, where where $\tilde{k}_{A \& B}$ is the rate of which a focal A particle is hit by B particles.


![](https://cdn.mathpix.com/cropped/2025_03_18_5e6b0bbe0ca706fbe34cg-4.jpg?height=669&width=1490&top_left_y=1777&top_left_x=259)


Since $x_{A}-x_{B}$ diffuses w/ twice the diffusion constant, we have that the number of collisions per unit time of B particles with a focal A particle is given by


$$
4 \pi(2 D)(2 R) C_{B}=\widetilde{k}_{A\leftarrow B}
$$

$$
\Rightarrow \frac{\partial C_{C}}{\partial t}=k_{\text {diff }} C_{A} \cdot C_{B} 
$$

where 

$$
\begin{gathered}k_{\text {diff }}=16 \pi D R=16 \pi R \frac{k_{B} T}{6 \pi n R}=\frac{8 k_{B} T}{3 \eta} \\ q\end{gathered}
$$

$$
\approx 7 \cdot 10^9  \frac{1}{M \cdot s}
$$

Notes:
- Most enzymes are 100-1000 times slower
$\Rightarrow$ they require more than just one collision to form a product.
- ... could be modeled by assuming a finite absorption rate.
- Relation to Michaelis-Menten kinetics: Rate of product formation =$k_{max}[E][S]/(K_M+[S])$.

Comparing with the above equation, we see that the enzyme efficience is bounded by $k_{diff}$
$$
\left[\frac{k_{\max }}{K_{\mu}}\right]=\frac{1}{\text { Concentration time }}
$$

