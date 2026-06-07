<div align="center">

<img width="200" height="200" alt="eagar-tsai-logo" src="https://github.com/user-attachments/assets/301d7569-b2c6-4e01-abaf-660d69689e90" />

# eagar-tsai

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
![Python](https://img.shields.io/badge/python-3.11%2B-blue)
![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20macos%20%7C%20windows-lightgrey)

[![Tests](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/test.yml/badge.svg)](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/test.yml)
[![Lint](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/lint.yml/badge.svg)](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/lint.yml)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19837225.svg)](https://doi.org/10.5281/zenodo.19837225)

`eagar-tsai` is a Python library implementing the Eagar–Tsai moving heat source model to estimate melt pool dimensions (length, width, depth) for a scanning laser over a semi-infinite solid. Temperature fields are computed via a 1D integral; melt pool dimensions are extracted from the liquidus isotherm. Built-in plotting covers 2-D temperature field heatmaps, 3-D temperature volume rendering with PyVista, and power-velocity printability maps.

<p>
  <a href="https://github.com/ArroyaveLab/eagar-tsai/issues/new?labels=bug">Report a Bug</a> |
  <a href="https://github.com/ArroyaveLab/eagar-tsai/issues/new?labels=enhancement">Request a Feature</a> |
  <a href="https://ArroyaveLab.github.io/eagar-tsai/">Documentation</a>
</p>

</div>

---

## Installation

Pre-built binary wheels are published to PyPI for Python 3.11, 3.12, and 3.13 on Linux (x86-64, i686), macOS (x86-64 and Apple Silicon), and Windows (AMD64). If a matching wheel exists for your platform, no C compiler is needed.

```sh
# Recommended — using uv
uv add eagar-tsai

# Alternative — using pip
pip install eagar-tsai
```

> [!NOTE]
> If no pre-built wheel matches your platform (for example, a non-standard Linux architecture or a Python version outside the supported range), the package falls back to building from source. In that case a C compiler is required: GCC or Clang on Linux/macOS; MSVC Build Tools or MinGW-w64 on Windows.

---

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

material = MaterialProperties(
        liquidus_temperature=1700.0,
        thermal_conductivity=30.0,
        density=7800.0,
        specific_heat=700.0
)

domain = SimulationDomain(
        x_length_um=1200.0,
        y_length_um=1200.0,
        z_depth_um=1000.0,
        spatial_resolution_um=1.0,
)

result = compute_single_point(beam, material, domain)

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
    data=df,
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

| Column              | Unit | Notes                                                                            |
|---------------------|------|----------------------------------------------------------------------------------|
| `melt_length`       | m    |                                                                                  |
| `melt_width`        | m    |                                                                                  |
| `melt_depth`        | m    |                                                                                  |
| `melt_length_um`    | µm   |                                                                                  |
| `melt_width_um`     | µm   |                                                                                  |
| `melt_depth_um`     | µm   |                                                                                  |
| `peak_temperature`  | K    |                                                                                  |
| `min_temperature`   | K    |                                                                                  |
| `temperature_field` | —    | `TemperatureField` object; `None` on failure. Pass `return_field=False` to omit. |

---

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

material  = MaterialProperties(
        liquidus_temperature=1700.0,
        thermal_conductivity=30.0,
        density=7800.0,
        specific_heat=700.0
)

domain = SimulationDomain(
        x_length_um=320.0,
        y_length_um=110.0,
        z_depth_um=60.0,
        spatial_resolution_um=5.0,
)

result = compute_single_point(
        beam=beam,
        material=material,
        domain=domain
)

print(result.length_um)                        # melt pool length in µm
print(result.temperature_field.T_xy.shape)     # (ny, nx) — surface plane in Kelvin
print(result.temperature_field.T_xz.shape)     # (nz, nx) — depth cross-section in Kelvin

fig = result.plot(output="temperature_field.png") # equivalently: result.temperature_field.plot(output="temperature_field.png")
```

The standalone convenience function skips constructing the `MeltPoolResult` object explicitly:

```python
from eagar_tsai.plot import plot_temperature_field

fig = plot_temperature_field(beam, material, domain, output="temperature_field.png")
```

---

## 3-D Temperature Volume

`compute_temperature_volume` evaluates the full `(nx, ny, nz)` temperature array on the auto-sized domain and returns a `TemperatureVolume`. The `y` axis covers only the half-domain by symmetry; `mirror_y=True` (the default) reconstructs the full melt pool for rendering and export.

The 3-D computation evaluates the integral at every grid point in the volume. For interactive use, coarser spatial resolution (e.g. `spatial_resolution_um=5.0`) gives a 25× speed-up over the 1 µm default with minimal loss of shape accuracy.

```python
from eagar_tsai import BeamParameters, MaterialProperties, SimulationDomain
from eagar_tsai import compute_temperature_volume

beam = BeamParameters(
    beam_diameter=80e-6,
    power=250.0,
    velocity=0.5,
    absorptivity=0.59,
)

material = MaterialProperties(
    liquidus_temperature=3455.0,
    thermal_conductivity=23.75,
    density=18038.9,
    specific_heat=251.6,
)

domain = SimulationDomain(
    x_length_um=800.0,
    y_length_um=400.0,
    z_depth_um=300.0,
    spatial_resolution_um=5.0,
)

volume = compute_temperature_volume(beam, material, domain, workers=-1)

print(volume.T_xyz.shape)              # (nx, ny, nz) in Kelvin
print(volume.result.length_um)         # melt pool length from the 2-D auto-sizing step

volume.export_vti("temperature_volume.vti")   # export for ParaView / VisIt

fig = volume.plot_3d()                 # returns a matplotlib Figure
plotter = volume.plot_3d(return_plotter=True)  # interactive PyVista window
```

The standalone convenience function skips constructing the `TemperatureVolume` object explicitly:

```python
from eagar_tsai.plot import plot_temperature_field_3d

# Returns a matplotlib Figure by default
fig = plot_temperature_field_3d(beam, material, domain, workers=-1)

# Save the rendered image and also export a .vti file in one call
fig = plot_temperature_field_3d(
    beam, material, domain,
    workers=-1,
    output="volume.png",
    output_vti="volume.vti",
)

# Open an interactive PyVista window
plotter = plot_temperature_field_3d(beam, material, domain, workers=-1, return_plotter=True)
```

To compute 3-D volumes for multiple parameter sets from a DataFrame, use `compute_temperature_volumes`:

```python
import pandas as pd
from eagar_tsai import compute_temperature_volumes

df = pd.DataFrame({
    "velocity_m_s":              [0.5, 1.0],
    "power_w":                   [200.0, 300.0],
    "beam_diameter_m":           [100e-6, 100e-6],
    "absorptivity":              [0.35, 0.35],
    "liquidus_temperature_k":    [1700.0, 1700.0],
    "thermal_conductivity_w_mk": [30.0, 30.0],
    "density_kg_m3":             [7800.0, 7800.0],
    "specific_heat_j_kgk":       [700.0, 700.0],
})

volumes = compute_temperature_volumes(df, workers=-1)
volumes[0].export_vti("row0_volume.vti")
```

---

## Printability Maps

`compute_printability_map` sweeps laser power and scan speed over a regular grid, runs the Eagar–Tsai model at every point, and classifies each point into one of four defect regimes (keyhole porosity, lack of fusion, balling, or defect-free) using the five criteria from Sheikh et al. (2023). Each grid point is dispatched as an independent parallel task, so workers stay fully utilized even when isolated points require iterative domain expansion.

```python
from eagar_tsai import (
    MaterialProperties,
    PrintabilityParameters,
    SimulationDomain,
    compute_printability_map,
)
from eagar_tsai.plot import plot_printability_map

material = MaterialProperties(
    liquidus_temperature=1700.0,
    thermal_conductivity=30.0,
    density=7800.0,
    specific_heat=700.0,
)

process = PrintabilityParameters(
    beam_diameter_m=80e-6,
    absorptivity=0.35,
    layer_thickness_m=40e-6,
    hatch_spacing_m=90e-6,
)

domain = SimulationDomain(
    x_length_um=1200.0,
    y_length_um=1200.0,
    z_depth_um=1000.0,
    spatial_resolution_um=5.0,  # coarser grid, ~25x faster than 1 µm
)

# Raw data: one row per grid point with defect classification
df = compute_printability_map(
    process,
    material,
    power_range=(50.0, 400.0),
    velocity_range=(0.1, 3.0),
    n_power=50,
    n_velocity=50,
    domain=domain,
    workers=-1,  # use all CPU cores
)

print(df["defect"].value_counts())

# Or render directly as a color-coded figure
fig = plot_printability_map(
    process,
    material,
    power_range=(50.0, 400.0),
    velocity_range=(0.1, 3.0),
    n_power=50,
    n_velocity=50,
    domain=domain,
    workers=-1,
    output="printability_map.png",
)
```

---

## References

- T. W. Eagar and N.-S. Tsai, "Temperature Fields Produced by Traveling Distributed Heat Sources," *Welding Journal (Research Supplement)*, December 1983, pp. 346-s–354-s.
- C integrand reformulation: Sasha Rubenchik, LLNL, 2015.

---

## License

This project is licensed under the GNU GPLv3 License. See the [LICENSE](https://github.com/ArroyaveLab/eagar-tsai/blob/main/LICENSE) file for details.

---

## Citation

We are currently preparing a preprint for publication. If you use `eagar-tsai` in your research, please cite the following:

> Sarıtürk, D., & Vela, B. (2026). eagar-tsai. Zenodo. https://doi.org/10.5281/zenodo.19837225

BibTeX:

```bibtex
@software{sariturk_2026_19837225,
  author    = {Sarıtürk, Doğuhan and Vela, Brent},
  title     = {eagar-tsai},
  year      = 2026,
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19837225},
  url       = {https://doi.org/10.5281/zenodo.19837225},
}
```
