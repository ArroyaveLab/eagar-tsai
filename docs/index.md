<div align="center" markdown>

# eagar-tsai
<img width="200" height="200" alt="eagar-tsai-logo" src="img/logo.svg" />

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
![Python](https://img.shields.io/badge/python-3.11%2B-blue)
![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20macos%20%7C%20windows-lightgrey)

[![Tests](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/test.yml/badge.svg)](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/test.yml)
[![Lint](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/lint.yml/badge.svg)](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/lint.yml)

`eagar-tsai` is a Python library implementing the Eagar–Tsai moving heat source model to estimate melt pool dimensions (length, width, depth) for a scanning laser over a semi-infinite solid. Temperature fields are computed via a 1D integral; melt pool dimensions are extracted from the liquidus isotherm.

<p>
  <a href="https://github.com/ArroyaveLab/eagar-tsai/issues/new?labels=bug">Report a Bug</a> |
  <a href="https://github.com/ArroyaveLab/eagar-tsai/issues/new?labels=enhancement">Request a Feature</a>
</p>

</div>

## Overview

The model computes temperature fields produced by a Gaussian laser beam moving over a semi-infinite solid. Melt pool dimensions are extracted from the liquidus isotherm.

**Reference:** T. W. Eagar and N.-S. Tsai, "Temperature Fields Produced by Traveling Distributed Heat Sources," *Welding Journal (Research Supplement)*, December 1983, pp. 346-s–354-s.

**Integrand reformulation:** Sasha Rubenchik, LLNL, 2015.

## Features

- **Fast integration** via a compiled C extension exposed as a `LowLevelCallable`, enabling QUADPACK to call the integrand at C speed with zero Python overhead per evaluation; pure-Python fallback when the extension is unavailable
- **Batch DataFrame API** with optional multiprocessing parallelism
- **Immutable dataclasses** for beam, material, and domain parameters
- **Iterative domain expansion** — automatically enlarges the simulation grid if the melt pool touches a boundary
- **Temperature field access** — `compute_single_point` returns a `MeltPoolResult` that always includes the full 2-D surface and depth temperature planes as an embedded `TemperatureField`; `result.plot()` produces a two-panel heatmap figure

## Installation

```bash
# Recommended — using uv
uv add eagar-tsai

# Alternative — using pip
pip install eagar-tsai
```

## Quick Start

For a single parameter set, use `compute_single_point` directly:

```python
from eagar_tsai import BeamParameters, MaterialProperties, SimulationDomain, compute_single_point

beam = BeamParameters(
        beam_diameter=100e-6,
        power=200.0,
        velocity=0.5,
        absorptivity=0.35
)
mat = MaterialProperties(
        liquidus_temperature=1700.0,
        thermal_conductivity=30.0,
        density=7800.0,
        specific_heat=700.0
)
dom = SimulationDomain(
        x_length_um=1200.0,
        y_length_um=1200.0,
        z_depth_um=1000.0,
        spatial_resolution_um=1.0,
)

result = compute_single_point(beam, mat, dom)
print(f"Length: {result.length_um:.1f} µm")
print(f"Width:  {result.width_um:.1f} µm")
print(f"Depth:  {result.depth_um:.1f} µm")
```

To run multiple parameter sets in parallel, use `compute_melt_pool` with a DataFrame. `workers`, `chunk_size`, and `output_dir` are optional:

```python
import pandas as pd
from eagar_tsai import compute_melt_pool

df = pd.DataFrame({
    "velocity_m_s":              [0.5],
    "power_w":                   [200.0],
    "beam_diameter_m":           [100e-6],
    "absorptivity":              [0.35],
    "liquidus_temperature_k":    [1700.0],
    "thermal_conductivity_w_mk": [30.0],
    "density_kg_m3":             [7800.0],
    "specific_heat_j_kgk":       [700.0],
})

result = compute_melt_pool(
    df,
    workers=4,                # parallel worker processes (default: None, serial)
    chunk_size=50,            # rows per worker chunk (default: 50)
    output_dir="CalcFiles/",  # write per-chunk CSVs (default: None)
)
print(result[["melt_length_um", "melt_width_um", "melt_depth_um"]])
```

### Required Input Columns

| Column                      | Unit     | Description                      |
|-----------------------------|----------|----------------------------------|
| `velocity_m_s`              | m/s      | Scan velocity                    |
| `power_w`                   | W        | Laser power                      |
| `beam_diameter_m`           | m        | Beam diameter (2σ)               |
| `absorptivity`              | —        | Absorptivity (0, 1]              |
| `liquidus_temperature_k`    | K        | Liquidus temperature             |
| `thermal_conductivity_w_mk` | W/(m·K)  | Thermal conductivity at liquidus |
| `density_kg_m3`             | kg/m³    | Density                          |
| `specific_heat_j_kgk`       | J/(kg·K) | Specific heat at liquidus        |

### Output Columns Added to the DataFrame

| Column             | Unit |
|--------------------|------|
| `melt_length`      | m    |
| `melt_width`       | m    |
| `melt_depth`       | m    |
| `melt_length_um`   | µm   |
| `melt_width_um`    | µm   |
| `melt_depth_um`    | µm   |
| `peak_temperature` | K    |
| `min_temperature`  | K    |

## Temperature Field Visualization

`compute_single_point` returns a `MeltPoolResult` that always includes the full `TemperatureField` and a built-in plot method:

```python
from eagar_tsai import BeamParameters, MaterialProperties, SimulationDomain, compute_single_point

beam = BeamParameters(
        beam_diameter=100e-6,
        power=200.0,
        velocity=0.5,
        absorptivity=0.35
)
mat  = MaterialProperties(
        liquidus_temperature=1700.0,
        thermal_conductivity=30.0,
        density=7800.0,
        specific_heat=700.0
)

result = compute_single_point(beam, mat)
print(result.length_um)                        # melt pool length in µm
print(result.temperature_field.T_xy.shape)     # (ny, nx) — surface plane in Kelvin
print(result.temperature_field.T_xz.shape)     # (nz, nx) — depth cross-section in Kelvin

fig = result.plot(output="temperature_field.png") # equivalently: result.temperature_field.plot(output="temperature_field.png")
```
