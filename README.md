<div align="center">

# eagar-tsai

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
![Python](https://img.shields.io/badge/python-3.11%2B-blue)
![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20macos%20%7C%20windows-lightgrey)

[![Tests](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/test.yml/badge.svg)](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/test.yml)
[![Lint](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/lint.yml/badge.svg)](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/lint.yml)

`eagar-tsai` is a Python library implementing the <strong>Eagar–Tsai moving heat source model</strong> to estimate melt pool dimensions (length, width, depth) for a scanning laser over a semi-infinite solid. Temperature fields are computed via a 1D integral; melt pool dimensions are extracted from the liquidus isotherm.

<p>
  <a href="https://github.com/ArroyaveLab/eagar-tsai/issues/new?labels=bug">Report a Bug</a> |
  <a href="https://github.com/ArroyaveLab/eagar-tsai/issues/new?labels=enhancement">Request a Feature</a> |
  <a href="https://ArroyaveLab.github.io/eagar-tsai/">Documentation</a>
</p>

</div>

---

## Installation

> [!NOTE]
> A C compiler is required: GCC or Clang on Linux/macOS; MSVC Build Tools or MinGW-w64 on Windows. The C extension is compiled automatically during installation.

```sh
# Recommended — using uv
uv add eagar-tsai

# Alternative — using pip
pip install eagar-tsai
```

---

## Quick Start

`compute_melt_pool` requires a DataFrame with the columns below; `workers`, `chunk_size`, and `output_dir` are optional:

```python
import pandas as pd
from eagar_tsai import compute_melt_pool

df = pd.DataFrame({
    "velocity_m_s":            [0.5],
    "power_w":                 [200.0],
    "beam_diameter_m":         [100e-6],
    "absorptivity":            [0.35],
    "liquidus_temperature_k":  [1700.0],
    "thermal_conductivity_w_mk": [30.0],
    "density_kg_m3":           [7800.0],
    "specific_heat_j_kgk":    [700.0],
})

result = compute_melt_pool(df)
print(result[["melt_length_um", "melt_width_um", "melt_depth_um"]])
```

Commonly overridden parameters:

```python
result = compute_melt_pool(
    df,
    workers=4,      # parallel worker processes (default: 1, serial)
    chunk_size=50,  # rows per worker chunk (default: 50)
    output_dir="CalcFiles/",  # write per-chunk CSVs (default: None)
)
```

### Required Input Columns

| Column | Unit | Description                      |
|--------|------|----------------------------------|
| `velocity_m_s` | m/s | Scan velocity                    |
| `power_w` | W | Laser power                      |
| `beam_diameter_m` | m | Beam diameter (2σ)               |
| `absorptivity` | — | Absorptivity [0, 1]              |
| `liquidus_temperature_k` | K | Liquidus temperature             |
| `thermal_conductivity_w_mk` | W/(m·K) | Thermal conductivity at liquidus |
| `density_kg_m3` | kg/m³ | Density                          |
| `specific_heat_j_kgk` | J/(kg·K) | Specific heat at liquidus        |

### Output Columns Added to the DataFrame

| Column | Unit |
|--------|------|
| `melt_length` | m |
| `melt_width` | m |
| `melt_depth` | m |
| `melt_length_um` | µm |
| `melt_width_um` | µm |
| `melt_depth_um` | µm |
| `peak_temperature` | K |
| `min_temperature` | K |

---

## Physics

**Reference:** T. W. Eagar and N.-S. Tsai, "Temperature Fields Produced by Traveling Distributed Heat Sources," *Welding Journal (Research Supplement)*, December 1983, pp. 346-s–354-s.

**C integrand reformulation:** Sasha Rubenchik, LLNL, 2015.

Thermal diffusivity:
```
alpha = k / (rho * cp)
```

Non-dimensional parameter:
```
p = alpha / (v * sigma)
```

Prefactor:
```
Ts = (A * P) / (pi * (k/alpha) * sqrt(pi * alpha * v * sigma^3))
```

Temperature field at (x, y, z):
```
T = T0 + Ts * integral_0^inf f(t, x, y, z, p) dt
```

Integrand:
```
f(t,x,y,z,p) = 1 / ((4pt+1) * sqrt(t))
               * exp(-z^2/(4t) - (y^2 + (x-t)^2)/(4pt+1))
```

Ambient temperature: `T0 = 300 K`.

### Assumptions

- Semi-infinite solid.
- Constant material properties at liquidus.
- No melt flow, vaporization, or latent heat effects.
- Gaussian heat source with constant absorptivity.
- Steady-state moving source.

---

## License

This project is licensed under the GNU GPLv3 License. See the [LICENSE](./LICENSE) file for details.
