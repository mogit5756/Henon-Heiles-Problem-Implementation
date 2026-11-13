# Hénon–Heiles Escape Dynamics

This project uses numerical simulation to study motion and escape in the two-dimensional Hénon–Heiles system. The [Jupyter notebook](Henon_Heiles_Python_Implementation.ipynb) contains the implementation, experiments, plots, and saved outputs. This README summarizes the physics and methods and maps the analysis to the notebook sections.

## Physical model

For unit mass, the dimensionless potential is

$$
V(x,y)=\frac12(x^2+y^2)+x^2y-\frac13y^3.
$$

The equations of motion follow from $\ddot{\mathbf r}=-\nabla V$:

$$
\ddot{x}=-x-2xy,\qquad \ddot{y}=-y-x^2+y^2.
$$

The conserved mechanical energy is $E=\frac12(v_x^2+v_y^2)+V(x,y)$. The three equivalent saddle points have critical energy $E_c=1/6$. Below this energy the central energetically allowed region is closed; above it, the three escape channels become accessible.

## Numerical method

Trajectories are integrated with **Velocity Verlet**, a second-order method for position-dependent acceleration:

$$
\mathbf r_{n+1}=\mathbf r_n+\mathbf v_n\Delta t+\frac12\mathbf a_n\Delta t^2,
$$
$$
\mathbf v_{n+1}=\mathbf v_n+\frac12(\mathbf a_n+\mathbf a_{n+1})\Delta t.
$$

Particles start at the origin with $v_{x0}=\sqrt{2E}\cos\theta$ and $v_{y0}=\sqrt{2E}\sin\theta$. The notebook checks energy conservation and timestep convergence; halving the timestep is expected to reduce the global error by about four in the second-order convergence regime.

## Escape definition

An escape is recorded when a particle crosses $R=2.5$ while moving outwards. The crossing angle is assigned to one of three approximate channel sectors (Upper, Lower-left, Lower-right); escapes outside these sectors are retained as Unclassified. A trajectory that does not cross $R$ before the final time is labelled **not escaped by $T$**, not permanently trapped.

## Experiments and notebook guide

| Notebook sections | Contents |
|---|---|
| 1–5 | Equations, escape geometry, Velocity Verlet, initial conditions, and plotting/error helpers |
| 6–7 | Energy conservation, timestep convergence, and symmetry-direction trajectories |
| 8 | Poincaré sections: phase-space crossings at selected energies |
| 9–11 | Full angular sweep, escape channels/times, and transition counting |
| 12–14 | Threefold symmetry, SciPy DOP853 comparison, and boundary timestep sensitivity |
| 15–17 | Coarse and refined energy scans |
| 18–20 | Fine basin scans, resolution comparison, and uncertainty fractions |

## Main numerical results

### Numerical validation

For the confined case $E=0.08$, $T=100$, the maximum absolute energy errors at $N=5{,}000$, $10{,}000$, and $20{,}000$ are $9.19\times10^{-6}$, $2.30\times10^{-6}$, and $5.74\times10^{-7}$. The error ratios (about 4.0004 and 4.0001) support second-order convergence. Rotating the upper-channel trajectory by 120° and 240° reproduces the symmetry-related trajectories with coordinate differences of order $10^{-14}$. For the representative $E=0.18$, $\theta=90^\circ$ escape, DOP853 gives $t_{esc}\approx5.48554$, compared with 5.49 for Velocity Verlet at $N=10{,}000$.

### Escape-basin investigation

At $E=0.18$, $T=100$, $N=20{,}000$, and 1° angular spacing, 228 of 360 sampled directions escape by $T$ (76 per channel); 132 do not escape within the observation time. The mean escape time among escaped trajectories is about 29.86.

In the 22°–23° scan at 0.01° spacing, the transition count changes from 37 at $N=20{,}000$ to 50 at $N=40{,}000$. The outcome changes for 31 of 100 sampled angles when the timestep is halved, so this local basin classification is not yet numerically converged.

### Energy parameter study

The refined scan uses $T=100$, $N=20{,}000$, and 5° spacing (72 launch directions per energy):

| Energy | Escaped / sampled | Fraction escaped by $T$ |
|---:|---:|---:|
| 0.167 | 6 / 72 | 0.083 |
| 0.170 | 18 / 72 | 0.250 |
| 0.175 | 24 / 72 | 0.333 |
| 0.180 | 54 / 72 | 0.750 |
| 0.200 | 72 / 72 | 1.000 |

These are finite-time fractions for the selected angular grid, not universal probabilities.

### Poincare sections

The section is $x=0$ with $v_x>0$, and the plotted coordinates are $(y,v_y)$. Fifteen initial conditions are sampled at each energy with $v_{y0}=0$; positive $v_{x0}$ is calculated from the fixed-energy condition. The initial crossing at $t=0$ is discarded. The saved sections contain 350, 347, 340, and 266 crossings at $E=1/12$, $1/8$, $1/6$, and 0.18. Six of the 15 trajectories at $E=0.18$ escape during the section run.

## Running the notebook

Install dependencies and launch Jupyter:

```bash
pip install numpy matplotlib scipy jupyter
jupyter lab
```

Open `Henon_Heiles_Python_Implementation.ipynb` and run the cells from top to bottom.

## Limitations

- Escape radius and channel sectors are numerical conventions; changing them can alter escape times and labels.
- Non-escape within the finite observation time does not establish permanent confinement.
- Boundary classifications are timestep-sensitive. DOP853 is stable for the selected tested angles, but it disagrees with some Velocity Verlet classifications. The local basin map is not demonstrated to be numerically converged.
- The uncertainty-fraction curves are exploratory: no uncertainty exponent or fractal dimension is claimed, and Wada behaviour has not been tested.
- Poincaré sections use finite ensembles and finite integration time, so they are qualitative phase-space diagnostics.

## Future work

Further work could include more systematic boundary convergence, uncertainty-exponent analysis after convergence is established, a higher-order symplectic integrator, and faster parameter scans.
