# ET Model (Eagar–Tsai)

[![CI](https://github.com/your-org/eagar-tsai/actions/workflows/ci.yml/badge.svg)](https://github.com/your-org/eagar-tsai/actions/workflows/ci.yml)

Python library implementing the **Eagar–Tsai moving heat source model** to estimate melt pool dimensions (length, width, depth) for a scanning laser over a semi-infinite solid. Temperature fields are computed via a 1D integral; melt pool dimensions are extracted from the liquidus isotherm.

**Reference:** T. W. Eagar and N.-S. Tsai, "Temperature Fields Produced by Traveling Distributed Heat Sources," *Welding Journal (Research Supplement)*, December 1983, pp. 346-s–354-s.

**C integrand reformulation:** Sasha Rubenchik, LLNL, 2015.

## Installation

```bash
# With uv (recommended)
uv add eagar-tsai

# With pip
pip install eagar-tsai
```

The C extension is compiled automatically. You need a C compiler (GCC/Clang on Linux/macOS; MSVC Build Tools or MinGW-w64 on Windows).

## Quick Start

```python
import pandas as pd
from eagar_tsai import compute_melt_pool

df = pd.DataFrame({
    "velocity_m_s": [0.5],
    "power_w": [200.0],
    "beam_diameter_m": [100e-6],
    "absorptivity": [0.35],
    "liquidus_temperature_k": [1700.0],
    "thermal_conductivity_w_mk": [30.0],
    "density_kg_m3": [7800.0],
    "specific_heat_j_kgk": [700.0],
})

result = compute_melt_pool(df)
print(result[["melt_length_um", "melt_width_um", "melt_depth_um"]])
```

## Required Input Columns

| Column | Unit | Description |
|--------|------|-------------|
| `velocity_m_s` | m/s | Scan velocity |
| `power_w` | W | Laser power |
| `beam_diameter_m` | m | Beam diameter (2σ) |
| `absorptivity` | — | Absorptivity (0, 1] |
| `liquidus_temperature_k` | K | Liquidus temperature |
| `thermal_conductivity_w_mk` | W/(m·K) | Thermal conductivity at liquidus |
| `density_kg_m3` | kg/m³ | Density |
| `specific_heat_j_kgk` | J/(kg·K) | Specific heat at liquidus |

## Output Columns

| Column | Unit |
|--------|------|
| `melt_length`, `melt_width`, `melt_depth` | m |
| `melt_length_um`, `melt_width_um`, `melt_depth_um` | µm |
| `peak_temperature`, `min_temperature` | K |

## Development Setup

```bash
# Clone and install in editable mode with all dev tools
git clone https://github.com/your-org/eagar-tsai
cd eagar-tsai
uv sync --extra dev

# Verify C extension compiled
python -c "from eagar_tsai._integrand_ext import get_integrand_capsule; print(get_integrand_capsule())"

# Run fast tests
uv run pytest tests/ -m "not slow"

# Run all tests
uv run pytest tests/

# Lint + format check
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/

# Type check
uv run ty check src/

# Install pre-commit hooks
uv run pre-commit install

# Docs preview
uv run --extra docs mkdocs serve
```

## Physics

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

## Assumptions

- Semi-infinite solid.
- Constant material properties at liquidus.
- No melt flow, vaporization, or latent heat effects.
- Gaussian heat source with constant absorptivity.
- Steady-state moving source.
