# Methodology

## Background

The Eagar–Tsai model (1983) provides an analytical solution for the steady-state temperature field produced by a Gaussian laser beam moving at constant velocity over a semi-infinite solid. It is widely used in additive manufacturing and laser welding research to estimate melt pool geometry without the cost of full finite-element simulations.

**References:**

- T. W. Eagar and N.-S. Tsai, "Temperature Fields Produced by Traveling Distributed Heat Sources," *Welding Journal (Research Supplement)*, December 1983, pp. 346-s–354-s.
- C integrand reformulation: Sasha Rubenchik, LLNL, 2015.

---

## Assumptions

- Semi-infinite solid (no boundaries other than the top surface).
- Constant material properties evaluated at the liquidus temperature.
- Gaussian heat source with constant absorptivity.
- Steady-state moving source (the temperature field moves with the beam).
- No melt flow, vaporization, or latent heat effects.

---

## Governing Equations

### Thermal diffusivity

Material thermal diffusivity is computed from the three input properties:

```
alpha = k / (rho * cp)
```

where `k` is thermal conductivity (W/(m·K)), `rho` is density (kg/m³), and `cp` is specific heat (J/(kg·K)).

### Non-dimensional parameter

The model uses a single non-dimensional parameter that captures the ratio of diffusive to advective transport:

```
p = alpha / (v * sigma)
```

where `v` is the scan velocity (m/s) and `sigma = sqrt(2) * (d / 2)` is the Gaussian width derived from beam diameter `d`.

### Temperature prefactor

The overall temperature scale is set by:

```
Ts = (A * P) / (pi * (k / alpha) * sqrt(pi * alpha * v * sigma^3))
```

where `A` is absorptivity and `P` is laser power (W).

### Temperature field

The temperature at any point `(x, y, z)` in the frame co-moving with the beam is:

```
T(x, y, z) = T0 + Ts * integral_0^inf  f(t, x, y, z, p)  dt
```

where `T0 = 298 K` is the ambient temperature and the integration variable `t` is a dimensionless time-like parameter.

### Integrand

```
f(t, x, y, z, p) =  1 / ((4*p*t + 1) * sqrt(t))
                  * exp(-z^2 / (4*t)  -  (y^2 + (x - t)^2) / (4*p*t + 1))
```

The integrand is evaluated numerically using `scipy.integrate.quad`. For performance, the integrand is implemented as a C extension (`_integrand_ext.c`) and passed to QUADPACK as a `LowLevelCallable`, eliminating Python overhead on every function evaluation.

---

## Melt Pool Extraction

The temperature field is evaluated on two planes:

- the **x–y plane** (z = 0, top surface) to obtain melt pool length and half-width,
- the **x–z plane** (y = 0, centerline) to obtain melt pool depth.

The melt pool boundary is the liquidus isotherm `T = T_liquidus`. The three dimensions are extracted as:

| Dimension | Definition |
|-----------|-----------|
| Length | Extent of `T >= T_liquidus` along x |
| Width | 2 × half-extent of `T >= T_liquidus` along y at the surface |
| Depth | Extent of `T >= T_liquidus` along z at the centerline |

If the melt pool reaches any domain boundary, the domain is automatically expanded and the computation is repeated (up to 20 iterations).
