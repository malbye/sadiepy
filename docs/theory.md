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
  The five correction notes at the end of the PDF fix Eqs. (6), (7), (14), (17) and
  (21); none affects Eqs. (8), (22) or (23), which are the ones used here.

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
  water. This is consistent with the paper's statement that the pure-water part
  deviates from IAPWS-95 by at most ±0.01 % (§2, p. 358): the fitted curve is
  flatter than real water near its maximum, but the density error stays of that
  order, far inside the stated ±0.1 % and irrelevant for wave damping. It does mean
  Eq. (1.1) should not be used to locate the density maximum.
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
and $x_0 \to \sqrt{y}$ (the shallow-water root of $x^2 = y$). Eq. (2.4) is the
explicit approximation of Eckart (1952), $L \approx L_0\sqrt{\tanh(2\pi h/L_0)}$ with
$L_0 = gT^2/2\pi$, rewritten in terms of $x$ and $y$ (named in the project spec §2.2;
`TODO(verify)`: primary source not yet in `refs/`). It is used only as a starting
point, so it does not affect the converged result. In between the limits it is
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

This is Zhang et al. (2023) Eq. (8), p. 3 (`refs/Zhang_et_al_2023.pdf`), which
states it without derivation.

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

These are Zhang et al. (2023) Eqs. (5)–(7), p. 3 (`refs/Zhang_et_al_2023.pdf`),
where the scales are written $u_s$ and $l_s$ and Ur is given in the form
$8\pi^3 (H/L)/(kh)^3$ with the linear wavelength $L = (gT^2/2\pi)\tanh kh$.
Zhang et al. use $u = \bar u$ (2.9, their §2.2.2) and
$\ell \in \{$DBH, $l_e$, $l_{ew}$, $r_v\}$ (their §2.2.3). The two forms of the Ursell number
are identical: with $kh = 2\pi h/L$,

$$
\frac{8\pi^3 (H/L)}{(kh)^3} = \frac{8\pi^3 H}{L}\cdot\frac{L^3}{8\pi^3 h^3} = \frac{H L^2}{h^3}.
$$

### 2.10 Validity checks

Linear theory assumes small steepness, $H/L \ll 1$, and small relative height,
$H/h \ll 1$. The code emits `ValidityWarning` (never stops silently) when:

| Check | Condition for warning | Source |
|---|---|---|
| Depth-limited breaking | $H/h > 0.78$ | Dean & Dalrymple (1991), Eq. (4.122a), p. 113 (McCowan, 1894), restated as Eq. (12.11), p. 335 for mild slopes; `refs/Dean_Dalrymple_1991.pdf` |
| Steepness-limited breaking | $H/L > 0.142\tanh kh$ | Miche (1944) form, chosen by the author (2026-11-23). `TODO(verify)`: citable source for the tanh form; see note below |
| Nonlinearity | $\mathrm{Ur} > \mathrm{Ur}_{\max}$ | `TODO(verify)`: default $\mathrm{Ur}_{\max}$ and its source; options to be presented at the M1 checkpoint |

$\mathrm{Ur}_{\max}$ is a documented, user-settable parameter.

**Depth-limited breaking.** Dean & Dalrymple (1991, p. 113) give McCowan's criterion
$H_b = \kappa h_b$ with $\kappa = 0.78$, and note on p. 335 that beach slopes steeper
than about 1:40 increase the ratio substantially. On the flat bed assumed here,
0.78 applies.

**Steepness-limited breaking (open).** Dean & Dalrymple (1991), Eq. (12.12),
p. 335, give only the *deep-water* limit $(H_b/L_o)_{\max} = 0.142$, with
$L_o = 1.2\,(gT^2/2\pi)$ "for breaking waves, including nonlinear effects". The
book does not state the depth-dependent form $H/L > 0.142\tanh kh$ (Miche, 1944)
used in the spec. In deep water the two differ: with the linear wavelength,
the spec's form triggers at $H/L = 0.142$, whereas the book's limit corresponds to
$H/L_{\text{linear}} = 0.142 \times 1.2 \approx 0.170$. The author chose the
Miche form (2026-11-23): it applies at every depth and is the more conservative
of the two in deep water. Its source remains `TODO(verify)` until a reference that
states the tanh form is in `refs/`.

---

## 3. Wave damping by vegetation

`sadiepy.core.solver` predicts how a regular wave decays across a patch of
vegetation described as a vertical stack of layers. The framework is that of
Dalrymple, Kirby & Hwang (1984) as used by Zhang et al. (2023), Eqs. (1)–(4),
pp. 2–3 (`refs/Zhang_et_al_2023.pdf`), generalised here to several layers.

**Assumptions.** Flat bed of constant depth $h$ (so $k$, $c_g$ are constant along
the patch); linear waves (§2); rigid vegetation elements (flexibility is M3); the
orbital velocity is that of the undisturbed wave (§2.7), i.e. the vegetation does
not modify the velocity profile; drag is the only dissipation (bed friction and
breaking are removed by bare-bed subtraction, M4); the wave propagates in $+x$
starting at $x = 0$ with height $H_0$.

### 3.1 Energy balance

With no reflection and a steady wave field, energy conservation across the patch is

$$
\frac{d}{dx}\left(E\,c_g\right) = -\varepsilon_D ,
\tag{3.1}
$$

(Zhang et al. 2023, Eq. (1), p. 2), where $E c_g$ is the energy flux (2.12) and
$\varepsilon_D$ is the time-averaged rate of energy dissipation per unit bed area
(W/m²).

### 3.2 Drag dissipation

Each element feels a drag force per unit height (Morison drag)
$f = \tfrac{1}{2}\rho\,C_D\,b(s)\,u\,|u|$, where $b(s)$ is the element's frontal
width at height $s$. The rate of work done by the fluid on it is $f\,u$. The
inertia force is proportional to $\partial u/\partial t$, which is $90°$ out of
phase with $u$, so its work averages to zero over a period and it is omitted.

A layer $i$ occupies $s_{b,i} \le s \le s_{t,i}$ and has $N_i$ elements per m² of
bed, each with frontal width $b_i(s)$ and drag coefficient $C_{D,i}$. Only the
submerged part dissipates energy, $s_{b,i} \le s \le s^*_{t,i}$ with
$s^*_{t,i} = \min(s_{t,i}, h)$; a layer with $s_{b,i} \ge h$ contributes nothing.
Summing over elements and layers,

$$
\varepsilon_D = \sum_i \Big\langle \int_{s_{b,i}}^{s^*_{t,i}}
\tfrac{1}{2}\rho\,C_{D,i}\,N_i\,b_i(s)\,|u|^3\,ds \Big\rangle
\tag{3.2}
$$

(cf. Zhang et al. 2023, Eq. (2), p. 3, where $A$ plays the role of $N\,b$).

### 3.3 Period average of $|\cos\theta|^3$

With $u = u_m(s)\cos\theta$ (2.8), $\langle|u|^3\rangle = u_m^3\,\langle|\cos\theta|^3\rangle$.
By symmetry the four quarter-periods contribute equally, so

$$
\langle|\cos\theta|^3\rangle
= \frac{1}{2\pi}\int_0^{2\pi}|\cos\theta|^3\,d\theta
= \frac{4}{2\pi}\int_0^{\pi/2}\cos^3\theta\,d\theta
= \frac{2}{\pi}\cdot\frac{2}{3}
= \frac{4}{3\pi},
\tag{3.3}
$$

using $\int_0^{\pi/2}\cos^3\theta\,d\theta = \int_0^1 (1 - v^2)\,dv = 2/3$
(substituting $v = \sin\theta$). Eq. (3.2) becomes

$$
\varepsilon_D = \sum_i \varepsilon_i,
\qquad
\varepsilon_i = \frac{2}{3\pi}\,\rho\,C_{D,i}\,N_i
\int_{s_{b,i}}^{s^*_{t,i}} b_i(s)\,u_m(s)^3\,ds .
\tag{3.4}
$$

### 3.4 Constant-width layers

Let $S(x) = \sinh^3 x + 3\sinh x$. Then $S'(x) = 3\sinh^2 x\cosh x + 3\cosh x
= 3\cosh x\,(\sinh^2 x + 1) = 3\cosh^3 x$, so

$$
\int_{s_b}^{s_t}\cosh^3(ks)\,ds = \frac{S(k s_t) - S(k s_b)}{3k}.
\tag{3.5}
$$

With $u_m(s) = \dfrac{H\omega}{2\sinh kh}\cosh ks$ (2.8) and constant $b_i$,
Eq. (3.4) gives

$$
\varepsilon_i = \frac{2}{3\pi}\,\rho\,C_{D,i}\,N_i\,b_i
\left(\frac{H\omega}{2\sinh kh}\right)^3
\frac{S(k s^*_{t,i}) - S(k s_{b,i})}{3k}
= A_i\,H^3 ,
\tag{3.6}
$$

where $A_i$ (units W/m⁵ = kg m⁻³ s⁻³) does not depend on $H$ when $C_{D,i}$ is
constant.

**Variable width.** When $b_i(s)$ is a function (e.g. the inverted-triangle
canopy), the integral in (3.4) is evaluated numerically
(`scipy.integrate.quad`). For constant $b_i$ the numerical and analytic (3.5)
results must agree (M1 test 5).

### 3.5 Hyperbolic decay for constant $C_D$

If every $C_{D,i}$ is constant, $\varepsilon_D = A H^3$ with $A = \sum_i A_i$.
With $E = \rho g H^2/8$ (2.11) and $c_g$ constant, (3.1) becomes

$$
\frac{\rho g c_g}{8}\,\frac{d(H^2)}{dx} = \frac{\rho g c_g}{4}\,H\frac{dH}{dx} = -A H^3
\quad\Longrightarrow\quad
\frac{d}{dx}\left(\frac{1}{H}\right) = -\frac{1}{H^2}\frac{dH}{dx} = \frac{4A}{\rho g c_g}.
$$

Integrating from $H(0) = H_0$:

$$
H(x) = \frac{H_0}{1 + \beta x},
\qquad
\beta = \frac{4 A H_0}{\rho g c_g},
\tag{3.7}
$$

(Zhang et al. 2023, Eq. (3), p. 3), with $\beta$ in m⁻¹. Note that $\beta$ depends
on $H_0$ but not on $\rho$, since $A \propto \rho$.

### 3.6 Reduction to Zhang et al. Eq. (4)

For one constant-width layer on $[0, h_{vs}]$ with $b = D$, combining (3.6) and
(3.7):

$$
\beta = \frac{4H_0}{\rho g c_g}\cdot\frac{2}{3\pi}\,\rho\,C_D N D\,
\frac{\omega^3}{8\sinh^3 kh}\cdot\frac{S(kh_{vs})}{3k}
= \frac{C_D D N H_0\,\omega^3\,S(kh_{vs})}{9\pi\,g\,c_g\,k\sinh^3 kh}.
$$

Substituting $c_g = \dfrac{\omega}{2k}\cdot\dfrac{\sinh 2kh + 2kh}{\sinh 2kh}$ (2.7),
then $\omega^2/g = k\tanh kh$ (2.2) and
$\sinh 2kh\,\tanh kh/\sinh^3 kh = 2/\sinh kh$:

$$
\beta = \frac{4}{9\pi}\,C_D D N H_0\,k\,
\frac{\sinh^3(kh_{vs}) + 3\sinh(kh_{vs})}{(\sinh 2kh + 2kh)\sinh kh},
\tag{3.8}
$$

which is Zhang et al. (2023), Eq. (4), p. 3. The layered solver must reproduce
(3.8) to $10^{-10}$ for a single uniform layer (M1 test 3).

**Note on $D$ versus $\bar D$.** Zhang et al. write Eq. (4) with $\bar D$, "the
depth averaged projected area", defined for the mangrove tree as
$\bar D = \frac{1}{h}\int A(z)\,dz$ and equal to the diameter $d$ for a cylinder.
For an emergent cylinder these agree, but for a submerged element ($h_{vs} < h$)
$\frac{1}{h}\int A\,dz = d\,h_{vs}/h$, which would count the submergence twice
(it is already in $S(kh_{vs})$). The layered solver avoids the issue by
integrating the actual width profile $b_i(s)$ (3.4); how $\bar D$ enters Zhang's
back-calculated $C_D$ values matters when their closures are applied (M2) and is
an open question for the author.

### 3.7 Additivity of $\beta$ across layers

From (3.7), $\beta = \frac{4H_0}{\rho g c_g}\sum_i A_i = \sum_i \beta_i$ with
$\beta_i = 4 A_i H_0 / (\rho g c_g)$: **damping factors add** when every layer has
a constant $C_D$, because the dissipation rates add (3.4).

*Layer splitting.* Splitting a layer at $s_m$ into two layers with the same
$C_D$, $N$ and $b$ leaves $\beta$ unchanged, because the bracket in (3.6)
telescopes: $[S(ks_m) - S(ks_b)] + [S(ks_t) - S(ks_m)] = S(ks_t) - S(ks_b)$
(M1 test 4).

*Zhang et al. Eq. (18).* Zhang et al. (2023), Eq. (18), p. 9, write
$\frac{1}{1+\beta x} = \sum_i \frac{1}{1+\beta_i x}$ for canopy, stem and root. This
is inconsistent with additive dissipation and fails even at $x = 0$, where the
left side is 1 and the right side is 3. Their own Eqs. (19)–(20) (p. 9) compute
each $\beta_i$ from Eq. (4) with the layer's bounds, which is the additive form
used here (spec erratum E3).

### 3.8 Drag coefficients that depend on the wave

When $C_{D,i}$ depends on the local wave (through Re, KC or Ur, all functions of
$H$), $\varepsilon_D$ is no longer proportional to $H^3$ and (3.7) does not hold.
Writing (3.1) as $\frac{\rho g c_g}{4} H \frac{dH}{dx} = -\varepsilon_D(H)$:

$$
\frac{dH}{dx} = -\frac{4\,\varepsilon_D(H)}{\rho\,g\,c_g\,H},
\qquad H(0) = H_0,
\tag{3.9}
$$

which is integrated along the patch with `scipy.integrate.solve_ivp`. The solver
reports $H(x)$ at the requested points, an effective $\beta$ from a least-squares
fit of (3.7) to $H(x)$ with its residual, each layer's share of the dissipation,
the range of Re, KC and Ur encountered, and any validity warnings.

**Checks.** With constant $C_D$, (3.9) must reproduce (3.7) to $10^{-6}$
(M1 test 6). Energy is conserved: integrating (3.1),
$E c_g\big|_{0} - E c_g\big|_{B} = \int_0^B \varepsilon_D\,dx$ over a patch of
length $B$, to a relative error below $10^{-6}$ (M1 test 7).

### 3.9 Bulk versus per-layer drag coefficients

Eq. (3.4) allows a different $C_{D,i}$ in each layer ("per-layer" $C_D$). Published
closures such as those in Zhang et al. (2023), Table 3, instead give one **bulk**
$C_D$ for the whole plant, defined with a particular frontal area and velocity
scale. Using a bulk closure means setting $C_{D,i} = C_D$ in every layer, with
$C_D$ evaluated from the closure's own definitions of Re, KC or Ur. Mathematically
this is (3.4) with equal coefficients; the solver records which mode was used.
Whether a bulk $C_D$ from one paper can be combined with this package's frontal
areas depends on how that paper defined its area (see the $\bar D$ note in §3.6);
the closures themselves are M2.

---

## 4. Froude scaling

A model at geometric scale $1:\lambda$ (prototype length / model length $= \lambda$)
reproduces the prototype's gravity-wave dynamics when the Froude number
$\mathrm{Fr} = u/\sqrt{g\,\ell}$ is the same in both, with the same $g$ and, as is
usual, the same fluid density and viscosity. Writing $\lambda_q$ for the
prototype/model ratio of a quantity $q$:

| Quantity | Relation | Scale factor |
|---|---|---|
| Length $\ell$ ($H$, $h$, $L$, $d$, $x$) | definition | $\lambda$ |
| Velocity $u$ | $\mathrm{Fr}$ equal $\Rightarrow u \propto \sqrt{\ell}$ | $\lambda^{1/2}$ |
| Time $T$ | $\ell / u$ | $\lambda^{1/2}$ |
| Wavenumber $k$ | 1/length | $\lambda^{-1}$ |
| Damping factor $\beta$ | 1/length, **only if $C_D$ is the same at both scales** | $\lambda^{-1}$ |
| Force | $\rho\,u^2\ell^2$ | $\lambda^{3}$ |
| $\mathrm{Re} = u\ell/\nu$ | same $\nu$ | $\lambda^{1/2}\cdot\lambda = \lambda^{3/2}$ |
| $\mathrm{KC} = uT/\ell$ | | $\lambda^{1/2}\cdot\lambda^{1/2}/\lambda = 1$ |
| $\mathrm{Ur} = HL^2/h^3$ | | $\lambda^{3}/\lambda^{3} = 1$ |
| Elastic modulus $E$ for equal Cauchy number | Cauchy number $\propto \rho u^2/E$ (definition to be verified, M3) | $\lambda$ |

The dispersion relation (2.2) is itself Froude-consistent: scaling $h$ by
$\lambda$ and $T$ by $\lambda^{1/2}$ leaves $y = \omega^2 h / g$ unchanged, so $kh$
is unchanged and $k$ scales as $\lambda^{-1}$. Hence KC and Ur are preserved and Re
is distorted by $\lambda^{3/2}$ (M1 test 8). Because $C_D$ generally depends on Re,
a model's $C_D$ (and hence $\beta\lambda$) need not equal the prototype's; this scale
effect is why closures fitted at one scale are re-evaluated at another (spec §2.5). These are the scalings listed in the
project spec §2.7; the Cauchy-number row is used only from M3 and its exact
definition is still to be verified against Luhar & Nepf (2016).
