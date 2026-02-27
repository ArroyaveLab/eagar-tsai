"""Shared pytest fixtures for the eagar-tsai test suite."""

from __future__ import annotations

import pytest

from eagar_tsai import BeamParameters, MaterialProperties, SimulationDomain


@pytest.fixture
def steel_beam() -> BeamParameters:
    """A representative steel laser beam configuration."""
    return BeamParameters(
        beam_diameter=100e-6,  # 100 um
        power=200.0,  # 200 W
        velocity=0.5,  # 0.5 m/s
        absorptivity=0.35,
    )


@pytest.fixture
def steel_material() -> MaterialProperties:
    """Approximate liquidus-temperature properties for 316L stainless steel."""
    return MaterialProperties(
        liquidus_temperature=1700.0,  # K
        thermal_conductivity=30.0,  # W/(m K)
        density=7800.0,  # kg/m^3
        specific_heat=700.0,  # J/(kg K)
    )


@pytest.fixture
def small_domain() -> SimulationDomain:
    """A coarse, small domain for fast tests (not physically accurate)."""
    return SimulationDomain(
        x_length_um=300.0,
        y_length_um=300.0,
        z_depth_um=200.0,
        spatial_resolution_um=10.0,
    )
