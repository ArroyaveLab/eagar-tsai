<div align="center" markdown>

# eagar-tsai
<img width="200" height="200" alt="eagar-tsai-logo" src="img/logo.svg" />

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
![Python](https://img.shields.io/badge/python-3.11%2B-blue)
![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20macos%20%7C%20windows-lightgrey)

[![Tests](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/test.yml/badge.svg)](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/test.yml)
[![Lint](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/lint.yml/badge.svg)](https://github.com/ArroyaveLab/eagar-tsai/actions/workflows/lint.yml)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19837225.svg)](https://doi.org/10.5281/zenodo.19837225)

`eagar-tsai` is a Python library implementing the Eagar–Tsai moving heat source model to estimate melt pool dimensions (length, width, depth) for a scanning laser over a semi-infinite solid. Temperature fields are computed via a 1D integral; melt pool dimensions are extracted from the liquidus isotherm. Built-in plotting covers temperature field heatmaps and power–velocity printability maps.

<p>
  <a href="https://github.com/ArroyaveLab/eagar-tsai/issues/new?labels=bug">Report a Bug</a> |
  <a href="https://github.com/ArroyaveLab/eagar-tsai/issues/new?labels=enhancement">Request a Feature</a>
</p>

</div>

---

## Overview

The model computes temperature fields produced by a Gaussian laser beam moving over a semi-infinite solid. Melt pool dimensions are extracted from the liquidus isotherm.

**Reference:** T. W. Eagar and N.-S. Tsai, "Temperature Fields Produced by Traveling Distributed Heat Sources," *Welding Journal (Research Supplement)*, December 1983, pp. 346-s–354-s.

**Integrand reformulation:** Sasha Rubenchik, LLNL, 2015.

---

## Features

- **Fast integration** via a compiled C extension exposed as a `LowLevelCallable`, enabling QUADPACK to call the integrand at C speed with zero Python overhead per evaluation; pure-Python fallback when the extension is unavailable
- **Batch DataFrame API** with optional multiprocessing parallelism
- **Printability maps**: sweep laser power x scan speed and classify every grid point into keyhole, lack of fusion, balling, or defect-free; each point is an independent parallel task for full CPU utilization
- **Immutable dataclasses** for beam, material, and domain parameters
- **Iterative domain expansion** — automatically enlarges the simulation grid if the melt pool touches a boundary
- **Temperature field access** — `compute_single_point` returns a `MeltPoolResult` that always includes the full 2-D surface and depth temperature planes as an embedded `TemperatureField`; `result.plot()` produces a two-panel heatmap figure
- **3D temperature volume** — `plot_temperature_field_3d` computes the full volumetric temperature distribution and renders an interactive or off-screen 3-D visualization via PyVista, with an optional liquidus isotherm contour surface and VTI export

---

## Installation

Install `eagar-tsai` with `uv` or `pip`. The C extension is compiled automatically during installation; no separate build step is needed.

```bash
# Recommended — using uv
uv add eagar-tsai

# Alternative — using pip
pip install eagar-tsai
```

---

## Quick Start

The example below computes melt pool dimensions for a single set of laser and material parameters. Pass `BeamParameters`, `MaterialProperties`, and `SimulationDomain` to `compute_single_point` to get length, width, and depth at the liquidus isotherm.

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
        spatial_resolution_um=1.0
)

result = compute_single_point(beam, material, domain)

print(f"Length: {result.length_um:.1f} µm")
print(f"Width:  {result.width_um:.1f} µm")
print(f"Depth:  {result.depth_um:.1f} µm")
```

See the [Usage](usage.md) page for batch DataFrame processing, temperature field visualization, printability maps, and parallel processing.
