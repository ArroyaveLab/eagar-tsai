"""Tests for eagar_tsai.plot."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest

from eagar_tsai import BeamParameters, MaterialProperties, SimulationDomain, compute_single_point
from eagar_tsai.plot import plot_temperature_field

_BEAM = BeamParameters(beam_diameter=100e-6, power=200.0, velocity=0.5, absorptivity=0.35)
_MATERIAL = MaterialProperties(
    liquidus_temperature=1700.0,
    thermal_conductivity=30.0,
    density=7800.0,
    specific_heat=700.0,
)
_DOMAIN = SimulationDomain(x_length_um=50.0, y_length_um=30.0, z_depth_um=20.0, spatial_resolution_um=10.0)


@pytest.mark.slow
def test_result_plot_returns_figure() -> None:
    """result.plot() returns a matplotlib Figure."""
    import matplotlib.figure

    result = compute_single_point(_BEAM, _MATERIAL, _DOMAIN)
    fig = result.plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_temperature_field_plot_returns_figure() -> None:
    """result.temperature_field.plot() returns a matplotlib Figure."""
    import matplotlib.figure

    result = compute_single_point(_BEAM, _MATERIAL, _DOMAIN)
    fig = result.temperature_field.plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_plot_temperature_field_returns_figure() -> None:
    """plot_temperature_field() returns a matplotlib Figure."""
    import matplotlib.figure

    fig = plot_temperature_field(_BEAM, _MATERIAL, _DOMAIN)
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_plot_saves_file(tmp_path) -> None:
    """result.plot(output=...) saves the figure to the specified path."""
    result = compute_single_point(_BEAM, _MATERIAL, _DOMAIN)
    out = tmp_path / "field.png"
    fig = result.plot(output=str(out))
    assert out.exists()
    plt.close(fig)


@pytest.mark.slow
def test_plot_annotate_false_returns_figure() -> None:
    """plot_temperature_field() with annotate=False still returns a Figure."""
    import matplotlib.figure

    fig = plot_temperature_field(_BEAM, _MATERIAL, _DOMAIN, annotate=False)
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)
