# Eagar–Tsai Model

Python library implementing the **Eagar–Tsai moving heat source model** for estimating melt pool dimensions (length, width, depth) during laser powder bed fusion and directed energy deposition processes.

## Overview

The model computes temperature fields produced by a Gaussian laser beam moving over a semi-infinite solid. Melt pool dimensions are extracted from the liquidus isotherm.

**Reference:** T. W. Eagar and N.-S. Tsai, "Temperature Fields Produced by Traveling Distributed Heat Sources," *Welding Journal (Research Supplement)*, December 1983, pp. 346-s–354-s.

**Integrand reformulation:** Sasha Rubenchik, LLNL, 2015.

## Features

- **Fast integration** via a compiled C extension exposed as a `LowLevelCallable`, enabling QUADPACK to call the integrand at C speed with zero Python overhead per evaluation; pure-Python fallback when the extension is unavailable
- **Batch DataFrame API** with optional multiprocessing parallelism
- **Immutable dataclasses** for beam, material, and domain parameters
- **Iterative domain expansion** — automatically enlarges the simulation grid if the melt pool touches a boundary

## Quick Install

```bash
pip install eagar-tsai
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv add eagar-tsai
```

## Quick Example

```python
import pandas as pd
from eagar_tsai import compute_melt_pool

df = pd.DataFrame({
    "velocity_m_s": [0.5, 1.0],
    "power_w": [200.0, 300.0],
    "beam_diameter_m": [100e-6, 100e-6],
    "absorptivity": [0.35, 0.35],
    "liquidus_temperature_k": [1700.0, 1700.0],
    "thermal_conductivity_w_mk": [30.0, 30.0],
    "density_kg_m3": [7800.0, 7800.0],
    "specific_heat_j_kgk": [700.0, 700.0],
})

result = compute_melt_pool(df, workers=2)
print(result[["melt_length_um", "melt_width_um", "melt_depth_um"]])
```
