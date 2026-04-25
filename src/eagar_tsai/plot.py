"""Matplotlib-based temperature field visualization for the Eagar-Tsai model.

Example:
    from eagar_tsai import BeamParameters, MaterialProperties, SimulationDomain
    from eagar_tsai.plot import plot_temperature_field

    beam = BeamParameters(beam_diameter=100e-6, power=200.0, velocity=0.5, absorptivity=0.35)
    mat = MaterialProperties(liquidus_temperature=1700.0, thermal_conductivity=30.0,
                              density=7800.0, specific_heat=700.0)
    fig = plot_temperature_field(beam, mat, output="temperature_field.png")
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

if TYPE_CHECKING:
    from pathlib import Path

    import matplotlib.figure

    from ._types import BeamParameters, MaterialProperties, SimulationDomain, TemperatureField

__all__ = ["plot_temperature_field"]


def plot_temperature_field(
    beam: BeamParameters,
    material: MaterialProperties,
    domain: SimulationDomain | None = None,
    *,
    output: str | Path | None = None,
    annotate: bool = True,
) -> matplotlib.figure.Figure:
    """Compute and render the x-y surface and x-z depth temperature planes.

    Internally calls ``compute_single_point`` and passes the temperature field
    to ``_render_temperature_panels``. Use this function when you only need the
    figure. If you also want access to the raw temperature arrays or melt pool
    summary, call ``compute_single_point`` directly and then call
    ``result.plot()`` or ``result.temperature_field.plot()``.

    Args:
        beam: Laser beam and process parameters.
        material: Material thermal properties.
        domain: Spatial domain. Defaults to ``SimulationDomain(1200, 1200, 1000, 1)`` um.
        output: File path to save the figure (e.g. ``"field.png"``). Supports
            any format recognized by ``matplotlib.figure.Figure.savefig``.
            When ``None`` the figure is returned without saving.
        annotate: When ``True`` (default), overlay width and depth annotations on the respective panels.

    Returns:
        A ``matplotlib.figure.Figure`` containing two panels: the x-y surface
        heatmap (top) and the x-z depth cross-section (bottom).
    """
    from . import compute_single_point

    result = compute_single_point(beam, material, domain)
    return _render_temperature_panels(result.temperature_field, output=output, annotate=annotate)


def _render_temperature_panels(
    field: TemperatureField,
    *,
    output: str | Path | None = None,
    annotate: bool = True,
) -> matplotlib.figure.Figure:
    """Render a two-panel temperature field figure from a ``TemperatureField``.

    Called by ``plot_temperature_field`` and by ``TemperatureField.plot()``.

    Args:
        field: Pre-computed temperature field.
        output: Optional file path to save the figure.
        annotate: When ``True``, overlay width and depth annotations derived
            from ``field.melt_width_m`` and ``field.melt_depth_m``.

    Returns:
        A ``matplotlib.figure.Figure``.
    """
    x_um = field.x_range_um
    y_um = field.y_range_um
    z_um = field.z_range_um  # negative values = depth below surface
    t_xy = field.T_xy
    t_xz = field.T_xz
    liquidus = field.liquidus_temperature_k

    # Mirror y for full (symmetric) melt pool view in panel a
    t_xy_mirror = np.vstack((np.flipud(t_xy[1:, :]), t_xy))
    y_full_um = np.concatenate((-y_um[:0:-1], y_um))

    # Flip to positive-down depth convention for panel b
    z_depth_um = -z_um[::-1]
    t_xz_display = np.flipud(t_xz)

    temp_norm = Normalize(vmin=300.0, vmax=max(float(t_xy.max()), float(t_xz.max())))
    cmap = plt.get_cmap("cividis")

    fig, (ax_a, ax_b) = plt.subplots(2, 1, sharex=True, figsize=(4.5, 4.0))

    # Panel a: x-y surface plane
    im_a = ax_a.imshow(
        t_xy_mirror,
        extent=[x_um[0], x_um[-1], y_full_um[0], y_full_um[-1]],
        origin="lower",
        aspect="auto",
        cmap=cmap,
        norm=temp_norm,
    )
    x_mesh_xy, y_mesh_xy = np.meshgrid(x_um, y_full_um)
    ax_a.contour(x_mesh_xy, y_mesh_xy, t_xy_mirror, levels=[liquidus], colors=["#d04f1a"], linewidths=1.0)
    div_a = make_axes_locatable(ax_a)
    cax_a = div_a.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im_a, cax=cax_a).set_label("T (K)")
    ax_a.set_ylabel("y (µm)")

    if annotate and field.melt_width_m > 0.0:
        w_um = field.melt_width_m * 1e6
        ax_a.annotate(
            "",
            xy=(0.0, w_um / 2.0),
            xytext=(0.0, -w_um / 2.0),
            arrowprops={"arrowstyle": "<->", "lw": 0.8, "color": "white"},
        )
        ax_a.text(
            10.0,
            0.0,
            f"W = {w_um:.0f} µm",
            color="white",
            rotation=90,
            va="center",
            ha="left",
            bbox={"boxstyle": "round,pad=0.15", "fc": (0, 0, 0, 0.28), "ec": "none"},
        )

    # Panel b: x-z depth cross-section
    im_b = ax_b.imshow(
        t_xz_display,
        extent=[x_um[0], x_um[-1], 0.0, z_depth_um[-1]],
        origin="lower",
        aspect="auto",
        cmap=cmap,
        norm=temp_norm,
    )
    x_mesh_xz, z_mesh_xz = np.meshgrid(x_um, z_depth_um)
    ax_b.contour(x_mesh_xz, z_mesh_xz, t_xz_display, levels=[liquidus], colors=["#d04f1a"], linewidths=1.0)
    div_b = make_axes_locatable(ax_b)
    cax_b = div_b.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im_b, cax=cax_b).set_label("T (K)")
    ax_b.invert_yaxis()
    ax_b.axhline(0.0, color="white", linewidth=0.55, alpha=0.8)
    ax_b.set_xlabel("x (µm)")
    ax_b.set_ylabel("Depth (µm)")

    if annotate and field.melt_depth_m > 0.0:
        d_um = field.melt_depth_m * 1e6
        ax_b.annotate(
            "",
            xy=(0.0, d_um),
            xytext=(0.0, 0.0),
            arrowprops={"arrowstyle": "<->", "lw": 0.8, "color": "white"},
        )
        ax_b.text(
            10.0,
            d_um / 2.0,
            f"D = {d_um:.0f} µm",
            color="white",
            rotation=90,
            va="center",
            ha="left",
            bbox={"boxstyle": "round,pad=0.15", "fc": (0, 0, 0, 0.28), "ec": "none"},
        )

    fig.tight_layout()

    if output is not None:
        fig.savefig(output, bbox_inches="tight")

    return fig
