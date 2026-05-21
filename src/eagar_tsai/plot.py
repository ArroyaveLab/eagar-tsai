"""Matplotlib-based temperature field visualization for the Eagar-Tsai model.

Example:
    from eagar_tsai import BeamParameters, MaterialProperties, SimulationDomain
    from eagar_tsai.plot import plot_temperature_field

    beam = BeamParameters(beam_diameter=100e-6, power=200.0, velocity=0.5, absorptivity=0.35)
    mat = MaterialProperties(liquidus_temperature=1700.0, thermal_conductivity=30.0, density=7800.0, specific_heat=700.0)
    fig = plot_temperature_field(beam, mat, output="temperature_field.png")
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "axes.linewidth": 0.6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.3,
        "ytick.major.size": 2.3,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

if TYPE_CHECKING:
    from pathlib import Path

    import matplotlib.figure

    from ._types import BeamParameters, MaterialProperties, PrintabilityParameters, SimulationDomain, TemperatureField

__all__ = ["plot_printability_map", "plot_temperature_field"]

_DEFECT_ORDER: list[str] = ["keyhole", "lack_of_fusion", "balling", "defect_free"]
_DEFECT_COLORS: list[str] = ["#4878cf", "#d65f5f", "#6acc65", "#f0f0f0"]
_DEFECT_DISPLAY: dict[str, str] = {
    "keyhole": "Keyhole",
    "lack_of_fusion": "Lack of Fusion",
    "balling": "Balling",
    "defect_free": "Defect-Free",
}


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

    temp_norm = Normalize(vmin=298.0, vmax=max(float(t_xy.max()), float(t_xz.max())))
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
    ax_a.set_ylabel("Width (µm)")
    ax_a.set_xlabel("Length (µm)")

    if annotate and field.melt_width_m > 0.0:
        w_um = field.melt_width_m * 1e6
        x_width_um = x_um[int(np.argmax(np.sum(t_xy_mirror >= liquidus, axis=0)))]
        ax_a.annotate(
            "",
            xy=(x_width_um, w_um / 2.0),
            xytext=(x_width_um, -w_um / 2.0),
            arrowprops={"arrowstyle": "<->", "lw": 0.8, "color": "white"},
        )
        ax_a.text(
            x_width_um + 10.0,
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
    ax_b.set_xlabel("Length (µm)")
    ax_b.set_ylabel("Depth (µm)")

    if annotate and field.melt_depth_m > 0.0:
        d_um = field.melt_depth_m * 1e6
        x_depth_um = x_um[int(np.argmax(np.sum(t_xz_display >= liquidus, axis=0)))]
        ax_b.annotate(
            "",
            xy=(x_depth_um, d_um),
            xytext=(x_depth_um, 0.0),
            arrowprops={"arrowstyle": "<->", "lw": 0.8, "color": "white"},
        )
        ax_b.text(
            x_depth_um + 10.0,
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


def plot_printability_map(
    params: PrintabilityParameters,
    material: MaterialProperties,
    *,
    power_range: tuple[float, float] = (40.0, 400.0),
    velocity_range: tuple[float, float] = (0.05, 3.0),
    n_power: int = 50,
    n_velocity: int = 50,
    keyhole_wdr_threshold: float = 2.5,
    domain: SimulationDomain | None = None,
    workers: int | None = None,
    output: str | Path | None = None,
    show_data_points: bool = False,
) -> matplotlib.figure.Figure:
    """Compute and render a printability map over a laser power * scan speed grid.

    Calls ``compute_printability_map`` internally and renders the resulting
    defect classification as a colour-coded map with laser power on the Y-axis
    and scan speed on the X-axis.

    Args:
        params: Fixed process parameters (beam diameter, absorptivity, layer
            thickness, hatch spacing).
        material: Material thermal properties.
        power_range: ``(min_power_W, max_power_W)`` for the grid. Defaults to ``(40.0, 400.0)``.
        velocity_range: ``(min_velocity_m_s, max_velocity_m_s)`` for the grid. Defaults to ``(0.05, 3.0)``.
        n_power: Number of laser power grid points. Defaults to ``50``.
        n_velocity: Number of scan speed grid points. Defaults to ``50``.
        keyhole_wdr_threshold: Width-to-depth ratio threshold for the KH1 keyhole
            criterion. Defaults to ``2.5``.
        domain: Simulation domain. For large grids a coarser domain
            (e.g. ``SimulationDomain(1200, 1200, 1000, 5)``) reduces compute time.
        workers: Worker processes for parallel computation. ``None`` or ``1``
            runs serially; ``-1`` uses all available cores.
        output: File path to save the figure. When ``None`` the figure is returned
            without saving.
        show_data_points: When ``True``, overlay a scatter marker at every
            computed grid point. Defaults to ``False``.

    Returns:
        A ``matplotlib.figure.Figure`` with a single axes showing the printability
        map coloured by defect regime.
    """
    from ._api import compute_printability_map

    df = compute_printability_map(
        params,
        material,
        power_range=power_range,
        velocity_range=velocity_range,
        n_power=n_power,
        n_velocity=n_velocity,
        keyhole_wdr_threshold=keyhole_wdr_threshold,
        domain=domain,
        workers=workers,
    )

    code_map = {label: i for i, label in enumerate(_DEFECT_ORDER)}
    codes = np.array([code_map[d] for d in df["defect"]], dtype=float).reshape(n_velocity, n_power)

    velocities = np.linspace(velocity_range[0], velocity_range[1], n_velocity)
    powers = np.linspace(power_range[0], power_range[1], n_power)

    from scipy.ndimage import gaussian_filter

    n_classes = len(_DEFECT_ORDER)
    one_hot = np.zeros((n_classes, n_velocity, n_power))
    for i in range(n_classes):
        one_hot[i] = (codes == i).astype(float)
    render = np.stack([gaussian_filter(one_hot[i], sigma=1.0) for i in range(n_classes)])

    fig, ax = plt.subplots(figsize=(3.504, 3.504 + 0.45), constrained_layout=True)

    _PLOT_ORDER = [1, 3, 2, 0]  # lack_of_fusion - defect_free - balling - keyhole
    for i in _PLOT_ORDER:
        if np.any(codes == i):
            others_max = np.max([render[j] for j in range(n_classes) if j != i], axis=0)
            margin = render[i] - others_max
            ax.contourf(velocities, powers, margin.T, levels=[0.0, 2.0], colors=[_DEFECT_COLORS[i]], alpha=0.8)

    if show_data_points:
        _color_map = {label: _DEFECT_COLORS[i] for i, label in enumerate(_DEFECT_ORDER)}
        ax.scatter(
            df["velocity_m_s"],
            df["power_w"],
            s=6,
            c=[_color_map[d] for d in df["defect"]],
            edgecolors="black",
            linewidths=0.4,
            zorder=3,
        )

    legend_handles = [
        Patch(facecolor=_DEFECT_COLORS[i], edgecolor="0.4", linewidth=0.6, label=_DEFECT_DISPLAY[label])
        for i, label in enumerate(_DEFECT_ORDER)
        if label in set(df["defect"])
    ]
    ax.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=len(legend_handles),
        frameon=False,
        fontsize=8,
        columnspacing=0.8,
        handlelength=1.0,
        handletextpad=0.4,
    )

    ax.set_xlabel("Scan Speed (m/s)")
    ax.set_ylabel("Laser Power (W)")

    if output is not None:
        fig.savefig(output, bbox_inches="tight")

    return fig
