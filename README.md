# README: Lee-Wave Energy Sinks in Bottom-Intensified Flow

This repository contains the numerical configurations required to reproduce the results presented in **Wu et al. (JPO, under revision)**. The study investigates the partition of lee-wave energy sinks in bottom-intensified flows, with particular emphasis on reabsorption into the mean flow versus explicit and indirect dissipation.

---

## Scientific Overview

This research identifies **Parametric Subharmonic Instability (PSI)** as a significant but previously overlooked nonlinear energy cascade that enhances lee-wave dissipation at the expense of reabsorption.

### Key Findings

* **PSI as a Major Sink:** Nonlinear energy transfer via PSI cascades energy from twice the Coriolis frequency toward the inertial frequency and high vertical wavenumbers, leading to indirect dissipation.

* **Energy Partition:** For typical broadband topographic height spectra at mid- and high latitudes, the net dissipative fraction is **90–95%**, while reabsorption accounts for only **5–10%**.

* **New Parameterization:** A new dissipative fraction parameterization, $y$, is proposed as a function of the normalized lee-wave generation frequency $\alpha\equiv|kU_0/f|$:

$$
y(\alpha) =
\begin{cases}
-0.1\alpha + 1.2, & \text{if } \alpha > 2 \\
1, & \text{if } \alpha \le 2
\end{cases}
$$

* **Negligible Radiation:** Free-wave radiation is found to be a negligible energy sink (~1%) and cannot explain the observed turbulence shortfall.

---

## Model Configurations

The study utilizes two independent numerical models to ensure that the findings are insensitive to vertical discretization, forcing, and dissipation schemes.

### 1. Process Study Ocean Model (PSOM)

* **Setup:** Configured with a bottom-intensified, laterally confined jet over sinusoidal topography.

* **Grid:** Terrain-following vertical discretization with 4 m resolution to align bottom cells with topography.

* **Viscosity:** Explores two extreme scenarios:
  * High viscosity: $4 \times 10^{-3} \ \mathrm{m}^2\mathrm{s}^{-1}$ to ensure stationarity.
  * Low viscosity: $10^{-4} \ \mathrm{m}^2\mathrm{s}^{-1}$ to permit wave–wave interactions.

### 2. MIT General Circulation Model (MITgcm)

* **Setup:** Extends different lee-wave generation regimes with topographic Froude numbers $Fr_t$ ranging from 0.05 (linear) to 0.5 (weakly nonlinear).

* **Grid:** Employs a partially filled-cell representation to handle steeper topography.

* **Mixing:** Utilizes the turbulence parameterization of Klymak et al. (2010).

---

## Repository Structure

* `/PSOM_leewaves` — Configuration files for lee waves over sinusoidal topography.

* `/MITgcm_test9_mono` — Sensitivity tests for higher $Fr_t$ cases.

---

## Contact

**Yue Cynthia Wu**  
Department of Naval Architecture and Marine Engineering  
University of Michigan, Ann Arbor  

Email: `ywuocean@umich.edu`
