# Getting Started

## Installation

### With uv (recommended)

```bash
uv add eagar-tsai
```

### With pip

```bash
pip install eagar-tsai
```

The C extension is compiled automatically during installation. You need a C compiler:

- **Linux/macOS**: GCC or Clang (usually pre-installed)
- **Windows**: MSVC Build Tools or MinGW-w64

## Input Data Format

`compute_melt_pool` accepts a `pandas.DataFrame` with the following columns:

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

The result DataFrame contains all input columns plus:

| Column | Unit | Description |
|--------|------|-------------|
| `melt_length` | m | Melt pool length |
| `melt_width` | m | Melt pool width (full, 2× half-width) |
| `melt_depth` | m | Melt pool depth |
| `melt_length_um` | µm | Melt pool length |
| `melt_width_um` | µm | Melt pool width |
| `melt_depth_um` | µm | Melt pool depth |
| `peak_temperature` | K | Peak temperature in domain |
| `min_temperature` | K | Minimum temperature in domain |

## Basic Usage

```python
import pandas as pd
from eagar_tsai import compute_melt_pool

df = pd.read_excel("my_parameters.xlsx")

# Rename columns to match required names if needed
df = df.rename(columns={
    "Velocity_m/s": "velocity_m_s",
    "Power": "power_w",
    "Beam_diameter_m": "beam_diameter_m",
    "Absorptivity": "absorptivity",
    "T_liquidus": "liquidus_temperature_k",
    "thermal_cond_liq": "thermal_conductivity_w_mk",
    "Density_kg/m3": "density_kg_m3",
    "Cp_J/kg": "specific_heat_j_kgk",
})

result = compute_melt_pool(df, workers=4, chunk_size=10)
result.to_csv("melt_pool_results.csv", index=False)
```

## Single-Point Computation

For a single parameter set, use `compute_single_point` directly. It returns a `MeltPoolResult` that always includes the full `TemperatureField`.

```python
from eagar_tsai import (
    BeamParameters, MaterialProperties, SimulationDomain, compute_single_point
)

beam = BeamParameters(
    beam_diameter=100e-6,         # m
    power=200.0,                  # W
    velocity=0.5,                 # m/s
    absorptivity=0.35,            # —
)
material = MaterialProperties(
    liquidus_temperature=1700.0,  # K
    thermal_conductivity=30.0,    # W/(m·K)
    density=7800.0,               # kg/m³
    specific_heat=700.0,          # J/(kg·K)
)
domain = SimulationDomain(
    x_length_um=1200.0,           # µm
    y_length_um=1200.0,           # µm
    z_depth_um=1000.0,            # µm
    spatial_resolution_um=1.0,    # µm
)

result = compute_single_point(beam, material, domain)
print(f"Length: {result.length_um:.1f} µm")
print(f"Width:  {result.width_um:.1f} µm")
print(f"Depth:  {result.depth_um:.1f} µm")
```

## Temperature Field Visualization

`compute_single_point` returns a `MeltPoolResult` that always includes a `TemperatureField` with direct access to the raw arrays and a built-in plot method.

```python
from eagar_tsai import (
    BeamParameters, MaterialProperties, SimulationDomain, compute_single_point
)

beam = BeamParameters(beam_diameter=100e-6, power=200.0, velocity=0.5, absorptivity=0.35)
mat  = MaterialProperties(liquidus_temperature=1700.0, thermal_conductivity=30.0,
                          density=7800.0, specific_heat=700.0)
dom  = SimulationDomain(x_length_um=320, y_length_um=110, z_depth_um=60,
                        spatial_resolution_um=5)

result = compute_single_point(beam, mat, dom)

# Access raw arrays via the embedded TemperatureField
print(result.temperature_field.T_xy.shape)   # (ny, nx) — surface plane in Kelvin
print(result.temperature_field.T_xz.shape)   # (nz, nx) — depth cross-section in Kelvin

# Render a two-panel figure (x-y surface heatmap + x-z depth heatmap)
fig = result.plot(output="temperature_field.png")
# equivalently: result.temperature_field.plot(output="temperature_field.png")
```

The standalone convenience function skips constructing the `MeltPoolResult` object explicitly:

```python
from eagar_tsai.plot import plot_temperature_field

fig = plot_temperature_field(beam, mat, dom, output="temperature_field.png")
```

## Parallel Processing

```python
# Use all available CPU cores
result = compute_melt_pool(df, workers=-1)

# Use 8 workers, process 50 rows per chunk
result = compute_melt_pool(df, workers=8, chunk_size=50)
```

!!! note
    Setting `workers=None` or `workers=1` runs serially in the current process.
    This is recommended for small datasets or debugging.

## Saving Intermediate Results

```python
from pathlib import Path

result = compute_melt_pool(
    df,
    workers=4,
    chunk_size=50,
    output_dir=Path("calc_files"),  # saves ET_0000.csv, ET_0001.csv, ...
)
```
