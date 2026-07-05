"""Eagar–Tsai moving heat source model for melt pool dimension estimation.

Quick-start example::

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

Reference: T. W. Eagar and N.-S. Tsai, "Temperature Fields Produced by
Traveling Distributed Heat Sources," Welding Journal (Research Supplement),
December 1983.

Reformulation: Sasha Rubenchik, LLNL, 2015.
"""

from ._api import compute_melt_pool, compute_printability_map, compute_temperature_volumes
from ._core import compute_single_point, compute_temperature_volume
from ._types import (
    BeamParameters,
    MaterialProperties,
    MeltPoolResult,
    PrintabilityParameters,
    SimulationDomain,
    TemperatureField,
    TemperatureVolume,
)
from .plot import plot_printability_map, plot_temperature_field, plot_temperature_field_3d

__version__ = "0.4.1"
__all__ = [
    "BeamParameters",
    "MaterialProperties",
    "MeltPoolResult",
    "PrintabilityParameters",
    "SimulationDomain",
    "TemperatureField",
    "TemperatureVolume",
    "compute_melt_pool",
    "compute_printability_map",
    "compute_single_point",
    "compute_temperature_volume",
    "compute_temperature_volumes",
    "plot_printability_map",
    "plot_temperature_field",
    "plot_temperature_field_3d",
]
