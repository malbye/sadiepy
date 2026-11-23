# SADIE theory

Derivations and source equations behind every physics function in `sadiepy`.
Equations are numbered `(section.n)`; code docstrings cite these numbers. Each
equation taken from a paper also names the paper's file in `refs/` (see
[refs/README.md](../refs/README.md)) and its equation number and page there.

All quantities are in SI units unless stated otherwise.

---

## 1. Water properties

`sadiepy.units.Water(temperature, salinity)` returns the density $\rho$ and
kinematic viscosity $\nu$ of fresh water or seawater at atmospheric pressure.
Both come from the correlations recommended by Sharqawy, Lienhard & Zubair (2010)
(`refs/Sharqawy_et_al_2010.pdf`).

### 1.1 Variables and scales

| Symbol | Meaning | Unit in the correlations |
|---|---|---|
| $t$ | temperature, ITS-90 scale | °C |
| $S$ | salinity, reference-composition (absolute) scale | kg/kg |
| $\rho_{sw}$ | seawater density | kg/m³ |
| $\mu_w$, $\mu_{sw}$ | dynamic viscosity of pure water and of seawater | kg/(m·s) = Pa·s |

The paper states its recommended correlations in ITS-90 temperature and
reference-composition salinity (Sharqawy et al. 2010, §1, p. 355). Note that
**$S$ enters the equations in kg/kg**, not g/kg: seawater at 35 g/kg is
$S = 0.035$. Salinity measured on the Practical Salinity Scale $S_P$ converts as
$S_R = 1.00472\,S_P$ (Sharqawy et al. 2010, Table 1, Eq. (3), p. 355); the
difference is below 0.5 % and is not applied automatically.

### 1.2 Density

Sharqawy et al. (2010), Table 2, Eq. (8), p. 357 (`refs/Sharqawy_et_al_2010.pdf`):

$$
\rho_{sw}(t, S) = \left(a_1 + a_2 t + a_3 t^2 + a_4 t^3 + a_5 t^4\right)
+ \left(b_1 S + b_2 S t + b_3 S t^2 + b_4 S t^3 + b_5 S^2 t^2\right)
\tag{1.1}
$$

| Coefficient | Value | Coefficient | Value |
|---|---|---|---|
| $a_1$ | $9.999 \times 10^{2}$ | $b_1$ | $8.020 \times 10^{2}$ |
| $a_2$ | $2.034 \times 10^{-2}$ | $b_2$ | $-2.001$ |
| $a_3$ | $-6.162 \times 10^{-3}$ | $b_3$ | $1.677 \times 10^{-2}$ |
| $a_4$ | $2.261 \times 10^{-5}$ | $b_4$ | $-3.060 \times 10^{-5}$ |
| $a_5$ | $-4.657 \times 10^{-8}$ | $b_5$ | $-1.613 \times 10^{-5}$ |

- **Validity:** $0 < t < 180$ °C, $0 < S < 0.16$ kg/kg, at $p = 0.1$ MPa.
- **Accuracy:** ±0.1 % (Table 2), against the data of Isdale & Morris and of
  Millero & Poisson (Fig. 2, p. 359).
- **Pure-water limit:** at $S = 0$, Eq. (1.1) reduces to the first bracket, which
  was fitted to IAPWS-95 pure-water density at 0.1 MPa with a maximum deviation of
  ±0.01 % (§2, p. 358).
- **Why this correlation:** it is the paper's recommended atmospheric-pressure
  density equation, it uses ITS-90 and absolute salinity directly (no scale
  conversion), and its range covers any laboratory flume or field condition.
  The printed errata (last page of the PDF) correct only Eq. (6), which is not used.

### 1.3 Dynamic viscosity

Sharqawy et al. (2010), Table 5, Eqs. (22)–(23), p. 364
(`refs/Sharqawy_et_al_2010.pdf`):

$$
\mu_{sw}(t, S) = \mu_w(t)\left(1 + A S + B S^2\right)
\tag{1.2}
$$

$$
A = 1.541 + 1.998 \times 10^{-2}\, t - 9.52 \times 10^{-5}\, t^2,
\qquad
B = 7.974 - 7.561 \times 10^{-2}\, t + 4.724 \times 10^{-4}\, t^2
\tag{1.3}
$$

with the pure-water viscosity, based on IAPWS 2008 data,

$$
\mu_w(t) = 4.2844 \times 10^{-5}
+ \left[\,0.157\,(t + 64.993)^2 - 91.296\,\right]^{-1}
\tag{1.4}
$$

- **Validity:** $0 < t < 180$ °C, $0 < S < 0.15$ kg/kg.
- **Accuracy:** ±1.5 % (Table 5); Fig. 7 (p. 365) shows the data scatter.
- **Why this correlation:** it is the paper's recommended equation, the only one in
  Table 5 written directly in ITS-90 and absolute salinity, and it reduces to
  the pure-water value (1.4) at $S = 0$.

### 1.4 Kinematic viscosity

By definition,

$$
\nu(t, S) = \frac{\mu_{sw}(t, S)}{\rho_{sw}(t, S)}
\tag{1.5}
$$

with $\mu_{sw}$ from (1.2) and $\rho_{sw}$ from (1.1). The combined validity is
the intersection of the two ranges: $0 < t < 180$ °C, $0 \le S \le 0.15$ kg/kg.

### 1.5 Hand-check values

Fresh water at $t = 20$ °C, $S = 0$:

- Eq. (1.1): $999.9 + 0.4068 - 2.4648 + 0.1809 - 0.0075 = 998.02$ kg/m³
  (IAPWS-95: 998.21 kg/m³; difference 0.02 %, within the stated ±0.1 %).
- Eq. (1.4): $0.157 \times 84.993^2 - 91.296 = 1042.84$, so
  $\mu_w = 4.2844 \times 10^{-5} + 1/1042.84 = 1.0018 \times 10^{-3}$ Pa·s
  (Swindells et al. reference value at 20 °C: $1.002 \times 10^{-3}$ Pa·s, quoted in
  Sharqawy et al. 2010, p. 364).
- Eq. (1.5): $\nu = 1.0018 \times 10^{-3} / 998.02 = 1.0038 \times 10^{-6}$ m²/s.

These values become unit tests for `units.Water`.

### 1.6 Implementation notes

- **Units at the API:** `Water(salinity=...)` takes salinity in **g/kg** (the unit
  people usually quote: fresh water 0, typical seawater 35) and divides by 1000
  before evaluating Eqs. (1.1)–(1.3).
- **Range checks:** the paper's ranges use strict inequalities, but $S = 0$ is the
  intended pure-water limit (§1.2, §1.3), so `Water` treats the bounds as
  inclusive: $0 \le t \le 180$ °C, $0 \le S \le 160$ g/kg for density and
  $0 \le S \le 150$ g/kg for viscosity. Inputs outside these ranges emit
  `ValidityWarning` and still return a value; they are never silently clipped.
  Negative salinity is unphysical and raises `ValueError`.
- **Density maximum:** the pure-water part of Eq. (1.1) peaks at $t \approx 1.67$ °C
  (where $a_2 + 2a_3 t + 3a_4 t^2 + 4a_5 t^3 = 0$), not at the ≈ 4 °C maximum of real
  water. The resulting error is about 0.002 % (999.884 vs 999.900 kg/m³ at 4 and
  0 °C), far inside the stated ±0.1 % accuracy and irrelevant for wave damping.
  It does mean Eq. (1.1) should not be used to locate the density maximum.
- **Pressure:** both correlations are for atmospheric pressure (0.1 MPa). Flume and
  nearshore depths of a few metres change density by about $5 \times 10^{-4}$ %
  per metre of depth, which is negligible here.

---

## 2. Linear wave theory

`sadiepy.core.waves` describes a regular (monochromatic) wave of height $H$ and
period $T$ on a flat bed of still-water depth $h$, to first order in wave
steepness (linear, or Airy, theory). Everything below is derived from the
velocity potential (2.1); no empirical coefficients enter except in the
validity checks of §2.10.

**Gravity.** $g = 9.80665$ m/s², the standard acceleration of gravity
(3rd CGPM, 1901). This replaces the rounded 9.81 m/s² of the original spec
(author's decision, 2026-11-22); the difference is 0.03 %.

### 2.1 Coordinates

Internally, height is measured **upward from the bed**: $s = z + h$, so
$0 \le s \le h$, with the bed at $s = 0$ and the still-water level at $s = h$.
Papers usually use $z$ from $-h$ (bed) to $0$ (still-water level); the conversion
$s = z + h$ is applied only where data enter or leave the package. The wave
propagates in $+x$, and the phase is $\theta = kx - \omega t$.

### 2.2 Notation

| Symbol | Meaning | Unit |
|---|---|---|
| $H$, $a = H/2$ | wave height, amplitude | m |
| $T$, $\omega = 2\pi/T$ | period, angular frequency | s, rad/s |
| $k$, $L = 2\pi/k$ | wavenumber, wavelength | rad/m, m |
| $h$ | still-water depth | m |
| $c$, $c_g$ | phase speed, group speed | m/s |
| $\eta$ | free-surface elevation above still-water level | m |

### 2.3 Velocity potential and free surface

For an incompressible, irrotational flow, the velocity is $\nabla\phi$ and
$\nabla^2\phi = 0$. With a no-flow bed ($\partial\phi/\partial s = 0$ at $s = 0$)
and the boundary conditions linearised about $s = h$, the progressive-wave
solution is

$$
\phi(x, s, t) = \frac{g a}{\omega}\,\frac{\cosh ks}{\cosh kh}\,\sin\theta,
\qquad
\eta(x, t) = a\cos\theta .
\tag{2.1}
$$

The cosh dependence satisfies the bed condition because
$\frac{d}{ds}\cosh ks = k\sinh ks = 0$ at $s = 0$. The amplitude of $\phi$ follows
from the linearised dynamic free-surface condition $\eta = -\frac{1}{g}\,\partial\phi/\partial t$
at $s = h$:
$-\frac{1}{g}\cdot\frac{g a}{\omega}\cdot(-\omega)\cos\theta = a\cos\theta$. ✓

### 2.4 Dispersion relation

The linearised kinematic condition $\partial\eta/\partial t = \partial\phi/\partial s$
at $s = h$, combined with the dynamic condition, gives
$\partial^2\phi/\partial t^2 + g\,\partial\phi/\partial s = 0$ at $s = h$.
Substituting (2.1):

$$
-\omega^2 \cosh kh + g k \sinh kh = 0
\quad\Longrightarrow\quad
\omega^2 = g k \tanh kh .
\tag{2.2}
$$

**Limits.** Deep water ($kh \to \infty$, $\tanh kh \to 1$): $k \to \omega^2/g$.
Shallow water ($kh \to 0$, $\tanh kh \approx kh$): $\omega^2 \approx g h k^2$, so
$k \to \omega/\sqrt{gh}$. These are M1 test 1.

### 2.5 Solving the dispersion relation

Write $x = kh$ and $y = \omega^2 h / g$ (known). Eq. (2.2) becomes

$$
f(x) = x \tanh x - y = 0,
\qquad
f'(x) = \tanh x + x\,\mathrm{sech}^2 x .
\tag{2.3}
$$

$f$ is strictly increasing for $x > 0$ (both terms of $f'$ are positive) with
$f(0) = -y < 0$ and $f \to \infty$, so there is exactly one positive root.

**Initial guess.**

$$
x_0 = \frac{y}{\sqrt{\tanh y}}
\tag{2.4}
$$

This is exact in both limits: for $y \gg 1$, $\tanh y \to 1$ and $x_0 \to y$
(the deep-water root, since $\tanh x \to 1$); for $y \ll 1$, $\tanh y \approx y$
and $x_0 \to \sqrt{y}$ (the shallow-water root of $x^2 = y$). In between it is
within 5 % of the root (checked numerically for $10^{-6} \le y \le 10^{4}$), so Newton's method

$$
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
\tag{2.5}
$$

converges to machine precision in at most 3 iterations over that range. The solver stops when
$|f(x_n)|/y < 10^{-14}$ and raises an error if that has not happened after a
fixed maximum number of iterations, rather than returning an unconverged value.
M1 test 1 checks the relative residual of (2.2) is below $10^{-12}$ over a wide
range of $y$.

### 2.6 Wavelength, phase speed, group speed

$$
L = \frac{2\pi}{k}, \qquad c = \frac{\omega}{k} = \frac{L}{T}.
\tag{2.6}
$$

The group speed is $c_g = d\omega/dk$. Differentiating (2.2),
$2\omega\,d\omega = g\left(\tanh kh + kh\,\mathrm{sech}^2 kh\right) dk$, so

$$
c_g = \frac{g\left(\tanh kh + kh\,\mathrm{sech}^2 kh\right)}{2\omega}
= \frac{c}{2}\left(1 + \frac{kh\,\mathrm{sech}^2 kh}{\tanh kh}\right)
= n\,c,
\qquad
n = \frac{1}{2}\left(1 + \frac{2kh}{\sinh 2kh}\right),
\tag{2.7}
$$

using $g/\omega = c/\tanh kh$ (from 2.2) and
$\mathrm{sech}^2 x / \tanh x = 1/(\sinh x \cosh x) = 2/\sinh 2x$.
Limits: $n \to 1/2$ in deep water, $n \to 1$ in shallow water.

### 2.7 Orbital velocity and excursion

The horizontal velocity is $u = \partial\phi/\partial x$. From (2.1),
$u = \frac{g a k}{\omega}\frac{\cosh ks}{\cosh kh}\cos\theta$, and by (2.2)
$gk/\omega = \omega/\tanh kh$, so $u = u_m(s)\cos\theta$ with amplitude

$$
u_m(s) = a\,\omega\,\frac{\cosh ks}{\sinh kh} = \frac{H}{2}\,\omega\,\frac{\cosh ks}{\sinh kh}.
\tag{2.8}
$$

**Depth-averaged amplitude.** Averaging (2.8) over $0 \le s \le h$:

$$
\bar u = \frac{1}{h}\int_0^h u_m(s)\,ds
= \frac{a\omega}{h\sinh kh}\cdot\frac{\sinh kh}{k}
= \frac{\omega a}{k h}.
\tag{2.9}
$$

This is Zhang et al. (2023) Eq. (8) (to be checked against
`refs/Zhang_et_al_2023.pdf` once added).

**Orbital excursion.** A particle's horizontal displacement is
$\int u\,dt = -(u_m/\omega)\sin\theta$, so its amplitude is

$$
A_w(s) = \frac{u_m(s)}{\omega}.
\tag{2.10}
$$

### 2.8 Energy and energy flux

**Energy density** (per unit surface area, averaged over a period
$\langle\cdot\rangle$). Potential energy relative to still water:
$\langle \tfrac{1}{2}\rho g \eta^2 \rangle = \tfrac{1}{2}\rho g a^2 \langle\cos^2\theta\rangle = \tfrac{1}{4}\rho g a^2$.
Kinetic energy, with $w = \partial\phi/\partial s = a\omega\frac{\sinh ks}{\sinh kh}\sin\theta$:

$$
\Big\langle \int_0^h \tfrac{1}{2}\rho\,(u^2 + w^2)\,ds \Big\rangle
= \frac{\rho a^2\omega^2}{4\sinh^2 kh}\int_0^h(\cosh^2 ks + \sinh^2 ks)\,ds
= \frac{\rho a^2\omega^2}{4\sinh^2 kh}\cdot\frac{\sinh 2kh}{2k}
= \frac{\rho a^2 \omega^2}{4k\tanh kh}
= \tfrac{1}{4}\rho g a^2 ,
$$

using $\cosh^2 + \sinh^2 = \cosh 2ks$, $\sinh 2kh = 2\sinh kh\cosh kh$ and (2.2).
Kinetic and potential energy are equal, and

$$
E = \tfrac{1}{2}\rho g a^2 = \frac{\rho g H^2}{8}.
\tag{2.11}
$$

**Energy flux** (per unit crest width): the rate of work by the dynamic
pressure $p = -\rho\,\partial\phi/\partial t = \rho g a\frac{\cosh ks}{\cosh kh}\cos\theta$
on the moving fluid, to leading order

$$
F = \Big\langle\int_0^h p\,u\,ds\Big\rangle
= \frac{\rho g a^2\omega}{2\sinh kh\cosh kh}\int_0^h\cosh^2 ks\,ds
= \frac{\rho g a^2\omega}{\sinh 2kh}\cdot\frac{2kh + \sinh 2kh}{4k}
= E\,\frac{c}{2}\left(1 + \frac{2kh}{\sinh 2kh}\right)
= E\,c_g ,
\tag{2.12}
$$

using $\int_0^h \cosh^2 ks\,ds = h/2 + \sinh(2kh)/(4k)$. Eq. (2.12) is the flux in
the energy balance of §3.

### 2.9 Dimensionless groups

With an explicitly chosen velocity scale $u$ and length scale $\ell$
(the code requires both as arguments and records the choice in its output):

$$
\mathrm{Re} = \frac{u\,\ell}{\nu}, \qquad
\mathrm{KC} = \frac{u\,T}{\ell}, \qquad
\mathrm{Ur} = \frac{H L^2}{h^3}.
\tag{2.13}
$$

These are Zhang et al. (2023) Eqs. (5)–(7) (to be checked against
`refs/Zhang_et_al_2023.pdf`). Zhang et al. use $u = \bar u$ (2.9) and
$\ell \in \{$DBH, $l_e$, $l_{ew}$, $r_v\}$. The two forms of the Ursell number
are identical: with $kh = 2\pi h/L$,

$$
\frac{8\pi^3 (H/L)}{(kh)^3} = \frac{8\pi^3 H}{L}\cdot\frac{L^3}{8\pi^3 h^3} = \frac{H L^2}{h^3}.
$$

### 2.10 Validity checks

Linear theory assumes small steepness, $H/L \ll 1$, and small relative height,
$H/h \ll 1$. The code emits `ValidityWarning` (never stops silently) when:

| Check | Condition for warning | Source |
|---|---|---|
| Depth-limited breaking | $H/h > 0.78$ | Project spec §2.2. `TODO(verify)`: primary source for 0.78 (usually attributed to McCowan, 1894) |
| Steepness-limited breaking | $H/L > 0.142\tanh kh$ | Project spec §2.2. `TODO(verify)`: primary source (usually attributed to Miche, 1944) |
| Nonlinearity | $\mathrm{Ur} > \mathrm{Ur}_{\max}$ | `TODO(verify)`: default $\mathrm{Ur}_{\max}$ and its source; options to be presented at the M1 checkpoint |

$\mathrm{Ur}_{\max}$ is a documented, user-settable parameter. Neither breaking
coefficient is in `refs/` yet; both are flagged for the author to supply or confirm.
