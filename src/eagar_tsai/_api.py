"""Public DataFrame API for batch Eagar-Tsai melt pool computation.

This module provides compute_melt_pool, the primary entry point for
processing multiple process parameter combinations stored in a
pandas DataFrame.
"""

from __future__ import annotations

import concurrent.futures
import logging
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ._core import compute_single_point
from ._types import BeamParameters, MaterialProperties, PrintabilityParameters, SimulationDomain

__all__ = ["compute_melt_pool", "compute_printability_map"]

_logger = logging.getLogger(__name__)

_REQUIRED_COLUMNS: list[str] = [
    "velocity_m_s",
    "power_w",
    "beam_diameter_m",
    "absorptivity",
    "liquidus_temperature_k",
    "thermal_conductivity_w_mk",
    "density_kg_m3",
    "specific_heat_j_kgk",
]


def _process_row(
    row: pd.Series,
    domain: SimulationDomain | None,
    return_field: bool = False,
) -> dict[str, Any]:
    """Process a single DataFrame row and return output columns as a dict.

    Args:
        row: A pandas Series with the required input columns.
        domain: Optional custom simulation domain.
        return_field: When ``True``, include ``temperature_field`` in the returned dict.

    Returns:
        Dictionary with keys:
        ``melt_length``, ``melt_width``, ``melt_depth``,
        ``melt_length_um``, ``melt_width_um``, ``melt_depth_um``,
        ``peak_temperature``, ``min_temperature``, and optionally
        ``temperature_field`` when ``return_field=True``.
    """
    beam = BeamParameters(
        beam_diameter=float(row["beam_diameter_m"]),
        power=float(row["power_w"]),
        velocity=float(row["velocity_m_s"]),
        absorptivity=float(row["absorptivity"]),
    )
    material = MaterialProperties(
        liquidus_temperature=float(row["liquidus_temperature_k"]),
        thermal_conductivity=float(row["thermal_conductivity_w_mk"]),
        density=float(row["density_kg_m3"]),
        specific_heat=float(row["specific_heat_j_kgk"]),
    )
    result = compute_single_point(beam, material, domain, full_field=return_field)
    out: dict[str, Any] = {
        "melt_length": result.length,
        "melt_width": result.width,
        "melt_depth": result.depth,
        "melt_length_um": result.length_um,
        "melt_width_um": result.width_um,
        "melt_depth_um": result.depth_um,
        "peak_temperature": result.peak_temperature,
        "min_temperature": result.min_temperature,
    }
    if return_field:
        out["temperature_field"] = result.temperature_field
    return out


def _process_chunk(
    params: tuple[int, pd.DataFrame, SimulationDomain | None, Path | None, bool],
) -> pd.DataFrame:
    """Process a chunk of rows and optionally save to CSV.

    Args:
        params: Tuple of ``(chunk_index, chunk_df, domain, output_dir, return_field)``.

    Returns:
        The chunk DataFrame with output columns appended.
    """
    chunk_idx, chunk, domain, output_dir, return_field = params
    chunk = chunk.reset_index(drop=True)

    output_records: list[dict[str, Any]] = []
    for _, row in chunk.iterrows():
        try:
            output_records.append(_process_row(row, domain, return_field))
        except Exception as exc:
            _logger.error("Error in chunk %d, row %s: %s", chunk_idx, row.name, exc)
            nan_row: dict[str, Any] = {
                "melt_length": float("nan"),
                "melt_width": float("nan"),
                "melt_depth": float("nan"),
                "melt_length_um": float("nan"),
                "melt_width_um": float("nan"),
                "melt_depth_um": float("nan"),
                "peak_temperature": float("nan"),
                "min_temperature": float("nan"),
            }
            if return_field:
                nan_row["temperature_field"] = None
            output_records.append(nan_row)

    out_df = pd.concat([chunk, pd.DataFrame(output_records, index=chunk.index)], axis=1)

    if output_dir is not None:
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        csv_path = out_path / f"ET_{chunk_idx:04d}.csv"
        out_df.to_csv(csv_path, index=False)
        _logger.debug("Wrote chunk %d to %s", chunk_idx, csv_path)

    return out_df


def _validate_columns(data: pd.DataFrame) -> None:
    """Raise ValueError if any required column is missing.

    Args:
        data: Input DataFrame to validate.

    Raises:
        ValueError: Lists every missing column in the error message.
    """
    missing = [col for col in _REQUIRED_COLUMNS if col not in data.columns]
    if missing:
        raise ValueError(f"Input DataFrame is missing required column(s): {missing}.\nRequired columns: {_REQUIRED_COLUMNS}")


def compute_melt_pool(
    data: pd.DataFrame,
    *,
    domain: SimulationDomain | None = None,
    chunk_size: int = 50,
    workers: int | None = None,
    output_dir: Path | str | None = None,
    return_field: bool = True,
) -> pd.DataFrame:
    """Compute melt pool dimensions for every row in a DataFrame.

    Parameters are processed in chunks; each chunk can be dispatched to a
    separate process for parallelism. Results are appended as new columns to
    a copy of ``data`` and returned.

    Args:
        data: Input DataFrame.  Must contain the columns listed in _REQUIRED_COLUMNS.
        domain: Custom simulation domain.  If None, the default 1200 x 1200 x 1000 um domain is used for every row.
        chunk_size: Number of rows per chunk.  Larger values reduce multiprocessing overhead at the cost of coarser progress.
            Defaults to 50.
        workers: Worker processes to use.  ``1`` or ``None`` runs serially; ``-1`` uses all available cores.
        output_dir: If provided, each processed chunk is saved as a CSV file under this directory before results are concatenated.
        return_field: When ``True``, a ``temperature_field`` column is added to the output DataFrame containing
            the ``TemperatureField`` for each row (``None`` for rows that failed).

    Returns:
        A new DataFrame identical to ``data`` plus the output columns:
        ``melt_length``, ``melt_width``, ``melt_depth`` (m),
        ``melt_length_um``, ``melt_width_um``, ``melt_depth_um`` (um),
        ``peak_temperature``, ``min_temperature`` (K), and optionally
        ``temperature_field`` when ``return_field=True``.

    Raises:
        TypeError: If ``data`` is not a pandas DataFrame.
        ValueError: If any required column is absent from ``data``.
        ValueError: If ``workers`` is not a positive integer, ``-1``, or ``None``.

    Examples:
        ```python
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
        result = compute_melt_pool(df, workers=1)
        ```
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError(f"data must be a pandas DataFrame, got {type(data).__name__!r}")
    _validate_columns(data)
    if workers is not None and workers != -1 and workers < 1:
        raise ValueError(f"workers must be a positive integer, -1 (all cores), or None, got {workers!r}")

    out_dir = Path(output_dir) if output_dir is not None else None

    chunks = [data.iloc[i : i + chunk_size] for i in range(0, len(data), chunk_size)]
    params = [(idx, chunk, domain, out_dir, return_field) for idx, chunk in enumerate(chunks)]

    if workers is None or workers == 1:
        _logger.info("Running serially (%d chunk(s)).", len(chunks))
        results = [_process_chunk(p) for p in params]
    else:
        max_workers = None if workers == -1 else workers
        _logger.info("Running with %s worker(s), %d chunk(s).", "all" if workers == -1 else workers, len(chunks))
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(_process_chunk, params))

    return pd.concat(results, ignore_index=True)


def _compute_printability_point(
    args: tuple[BeamParameters, MaterialProperties, SimulationDomain],
) -> tuple[float, float, float]:
    """Compute melt pool dimensions for a single printability map grid point.

    Args:
        args: Tuple of (BeamParameters, MaterialProperties, SimulationDomain).

    Returns:
        Tuple of (length_um, width_um, depth_um). Returns NaN values on error.
    """
    beam, material, domain = args
    try:
        result = compute_single_point(beam, material, domain, full_field=False)
        return result.length_um, result.width_um, result.depth_um
    except Exception as exc:
        _logger.error("Error computing printability point (P=%.1f W, v=%.3f m/s): %s", beam.power, beam.velocity, exc)
        return float("nan"), float("nan"), float("nan")


def _classify_defect(
    length_um: float,
    width_um: float,
    depth_um: float,
    layer_thickness_um: float,
    hatch_spacing_um: float,
    keyhole_wdr_threshold: float = 2.5,
) -> tuple[str, bool, bool, bool, bool, bool]:
    """Classify a single melt pool point into a defect regime.

    Applies five physics-based criteria from Sheikh et al. (2023) in priority
    order: keyhole > lack of fusion > balling > defect-free.

    Args:
        length_um: Melt pool length in µm.
        width_um: Melt pool full width in µm.
        depth_um: Melt pool depth in µm.
        layer_thickness_um: Powder layer thickness in µm (LOF criteria).
        hatch_spacing_um: Hatch spacing in µm (LOF2 criterion).
        keyhole_wdr_threshold: Width-to-depth ratio below which keyholing is
            predicted (KH1). Default 2.5 per Sheikh et al.

    Returns:
        A tuple ``(defect, lof1, lof2, ball1, ball2, kh1)`` where ``defect``
        is one of ``"defect_free"``, ``"keyhole"``, ``"lack_of_fusion"``, or
        ``"balling"`` and the remaining five elements are the individual
        criterion flags.
    """
    if depth_um <= 0.0 or width_um <= 0.0:
        return "lack_of_fusion", True, True, False, False, False

    lof1 = depth_um <= layer_thickness_um
    lof2 = (hatch_spacing_um / width_um) ** 2 + layer_thickness_um / (layer_thickness_um + depth_um) >= 1.0
    ball1 = length_um / width_um >= 2.3
    ball2 = math.pi * width_um / length_um < math.sqrt(2.0 / 3.0)
    kh1 = width_um / depth_um <= keyhole_wdr_threshold

    if kh1:
        defect = "keyhole"
    elif lof1 or lof2:
        defect = "lack_of_fusion"
    elif ball1 or ball2:
        defect = "balling"
    else:
        defect = "defect_free"

    return defect, lof1, lof2, ball1, ball2, kh1


def compute_printability_map(
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
) -> pd.DataFrame:
    """Compute a printability map over a laser power * scan speed grid.

    Runs the Eagar-Tsai model at every (power, velocity) grid point and
    classifies each point using five defect criteria from Sheikh et al. (2023):
    LOF1, LOF2, Ball1, Ball2, and KH1. Points are labeled as one of
    ``"defect_free"``, ``"keyhole"``, ``"lack_of_fusion"``, or ``"balling"`` in
    priority order: keyhole > lack of fusion > balling > defect-free.

    Args:
        params: Fixed process parameters (beam diameter, absorptivity, layer
            thickness, hatch spacing).
        material: Material thermal properties.
        power_range: ``(min_power_W, max_power_W)`` for the grid. Defaults to ``(40.0, 400.0)``.
        velocity_range: ``(min_velocity_m_s, max_velocity_m_s)`` for the grid. Defaults to ``(0.05, 3.0)``.
        n_power: Number of power grid points. Defaults to ``50``.
        n_velocity: Number of velocity grid points. Defaults to ``50``.
        keyhole_wdr_threshold: Width-to-depth ratio threshold for KH1 keyhole
            criterion. Defaults to ``2.5``.
        domain: Simulation domain. Defaults to ``SimulationDomain(1200, 1200, 1000, 5)``
            (5 µm resolution), which is ~25x faster than the 1 µm default used
            by ``compute_melt_pool`` with negligible classification accuracy loss.
            Pass an explicit domain to override.
        workers: Worker processes for parallel computation. ``None`` or ``1``
            runs serially; ``-1`` uses all available cores. Each grid point is
            dispatched as an independent task, so workers stay fully utilized
            even when isolated points require iterative domain expansion.

    Returns:
        A DataFrame with one row per grid point containing:
        ``power_w``, ``velocity_m_s``, ``melt_length_um``, ``melt_width_um``,
        ``melt_depth_um``, ``defect`` (str), ``lof1``, ``lof2``, ``ball1``,
        ``ball2``, ``kh1`` (bool).

    Examples:
        ```python
        from eagar_tsai import MaterialProperties, PrintabilityParameters, compute_printability_map

        mat = MaterialProperties(
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
        df = compute_printability_map(process, mat, n_power=30, n_velocity=30, workers=-1)
        print(df["defect"].value_counts())
        ```
    """
    powers = np.linspace(power_range[0], power_range[1], n_power)
    velocities = np.linspace(velocity_range[0], velocity_range[1], n_velocity)
    pv_grid, vv_grid = np.meshgrid(powers, velocities)

    effective_domain = (
        domain
        if domain is not None
        else SimulationDomain(
            x_length_um=1200.0,
            y_length_um=1200.0,
            z_depth_um=1000.0,
            spatial_resolution_um=5.0,
        )
    )

    job_args: list[tuple[BeamParameters, MaterialProperties, SimulationDomain]] = [
        (
            BeamParameters(
                beam_diameter=params.beam_diameter_m,
                power=float(p),
                velocity=float(v),
                absorptivity=params.absorptivity,
            ),
            material,
            effective_domain,
        )
        for p, v in zip(pv_grid.ravel(), vv_grid.ravel(), strict=False)
    ]

    if workers is None or workers == 1:
        melt_results = [_compute_printability_point(a) for a in job_args]
    else:
        max_workers = None if workers == -1 else workers
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            melt_results = list(executor.map(_compute_printability_point, job_args))

    lengths_um = [r[0] for r in melt_results]
    widths_um = [r[1] for r in melt_results]
    depths_um = [r[2] for r in melt_results]

    layer_um = params.layer_thickness_um
    hatch_um = params.hatch_spacing_um

    classifications = [
        _classify_defect(ln, w, d, layer_um, hatch_um, keyhole_wdr_threshold)
        for ln, w, d in zip(lengths_um, widths_um, depths_um, strict=True)
    ]

    defects, lof1s, lof2s, ball1s, ball2s, kh1s = zip(*classifications, strict=True)

    return pd.DataFrame(
        {
            "power_w": pv_grid.ravel(),
            "velocity_m_s": vv_grid.ravel(),
            "melt_length_um": lengths_um,
            "melt_width_um": widths_um,
            "melt_depth_um": depths_um,
            "defect": list(defects),
            "lof1": list(lof1s),
            "lof2": list(lof2s),
            "ball1": list(ball1s),
            "ball2": list(ball2s),
            "kh1": list(kh1s),
        }
    )
