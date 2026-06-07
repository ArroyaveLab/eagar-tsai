"""Tests for eagar_tsai.plot."""

from __future__ import annotations

import matplotlib.figure
import matplotlib.pyplot as plt
import pytest

from eagar_tsai import compute_single_point, compute_temperature_volume
from eagar_tsai.plot import plot_printability_map, plot_temperature_field, plot_temperature_field_3d


@pytest.mark.slow
def test_result_plot_returns_figure(steel_beam, steel_material, tiny_domain) -> None:
    """result.plot() returns a matplotlib Figure."""
    import matplotlib.figure

    result = compute_single_point(steel_beam, steel_material, tiny_domain)
    fig = result.plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_temperature_field_plot_returns_figure(steel_beam, steel_material, tiny_domain) -> None:
    """result.temperature_field.plot() returns a matplotlib Figure."""
    import matplotlib.figure

    result = compute_single_point(steel_beam, steel_material, tiny_domain)
    fig = result.temperature_field.plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_plot_temperature_field_returns_figure(steel_beam, steel_material, tiny_domain) -> None:
    """plot_temperature_field() returns a matplotlib Figure."""
    import matplotlib.figure

    fig = plot_temperature_field(steel_beam, steel_material, tiny_domain)
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_plot_saves_file(tmp_path, steel_beam, steel_material, tiny_domain) -> None:
    """result.plot(output=...) saves the figure to the specified path."""
    result = compute_single_point(steel_beam, steel_material, tiny_domain)
    out = tmp_path / "field.png"
    fig = result.plot(output=str(out))
    assert out.exists()
    plt.close(fig)


@pytest.mark.slow
def test_plot_annotate_false_returns_figure(steel_beam, steel_material, tiny_domain) -> None:
    """plot_temperature_field() with annotate=False still returns a Figure."""
    import matplotlib.figure

    fig = plot_temperature_field(steel_beam, steel_material, tiny_domain, annotate=False)
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_plot_3d_returns_figure(volume_domain, steel_beam, steel_material) -> None:
    """TemperatureVolume.plot_3d() returns a matplotlib Figure by default."""
    vol = compute_temperature_volume(steel_beam, steel_material, volume_domain)
    fig = vol.plot_3d()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_plot_3d_saves_file(tmp_path, volume_domain, steel_beam, steel_material) -> None:
    """TemperatureVolume.plot_3d(output=...) saves the figure to disk."""
    vol = compute_temperature_volume(steel_beam, steel_material, volume_domain)
    out = tmp_path / "volume.png"
    fig = vol.plot_3d(output=str(out))
    assert out.exists()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_plot_temperature_field_3d_returns_figure(volume_domain, steel_beam, steel_material) -> None:
    """plot_temperature_field_3d() returns a matplotlib Figure by default."""
    fig = plot_temperature_field_3d(steel_beam, steel_material, volume_domain)
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_export_vti_creates_file(tmp_path, volume_domain, steel_beam, steel_material) -> None:
    """TemperatureVolume.export_vti() writes a .vti file to the given path."""
    vol = compute_temperature_volume(steel_beam, steel_material, volume_domain)
    out = tmp_path / "volume.vti"
    returned_path = vol.export_vti(str(out))
    assert out.exists()
    assert returned_path == out.resolve()


@pytest.mark.slow
def test_plot_printability_map_returns_figure(steel_material, printability_domain, printability_params) -> None:
    """plot_printability_map() returns a matplotlib Figure."""
    import matplotlib.figure

    fig = plot_printability_map(
        printability_params,
        steel_material,
        power_range=(100.0, 300.0),
        velocity_range=(0.5, 1.5),
        n_power=3,
        n_velocity=3,
        domain=printability_domain,
    )
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


@pytest.mark.slow
def test_plot_printability_map_saves_file(tmp_path, steel_material, printability_domain, printability_params) -> None:
    """plot_printability_map(output=...) saves the figure to the specified path."""
    out = tmp_path / "pmap.png"
    fig = plot_printability_map(
        printability_params,
        steel_material,
        power_range=(100.0, 300.0),
        velocity_range=(0.5, 1.5),
        n_power=3,
        n_velocity=3,
        domain=printability_domain,
        output=str(out),
    )
    assert out.exists()
    plt.close(fig)
