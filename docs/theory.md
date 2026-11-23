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
