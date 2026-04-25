"""Public DataFrame API for batch Eagar-Tsai melt pool computation.

This module provides compute_melt_pool, the primary entry point for
processing multiple process parameter combinations stored in a
pandas DataFrame.
"""

from __future__ import annotations

import concurrent.futures
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from ._core import compute_single_point
from ._types import BeamParameters, MaterialProperties

if TYPE_CHECKING:
    from ._types import SimulationDomain

__all__ = ["compute_melt_pool"]

_logger = logging.getLogger(__name__)

REQUIRED_COLUMNS: list[str] = [
    "velocity_m_s",
    "power_w",
    "beam_diameter_m",
    "absorptivity",
    "liquidus_temperature_k",
    "thermal_conductivity_w_mk",
    "density_kg_m3",
    "specific_heat_j_kgk",
]

_NAN_RESULT: dict[str, Any] = {
    "melt_length": float("nan"),
    "melt_width": float("nan"),
    "melt_depth": float("nan"),
    "melt_length_um": float("nan"),
    "melt_width_um": float("nan"),
    "melt_depth_um": float("nan"),
    "peak_temperature": float("nan"),
    "min_temperature": float("nan"),
}


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
    result = compute_single_point(beam, material, domain)
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
            nan_row: dict[str, Any] = _NAN_RESULT.copy()
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
    missing = [col for col in REQUIRED_COLUMNS if col not in data.columns]
    if missing:
        raise ValueError(f"Input DataFrame is missing required column(s): {missing}.\nRequired columns: {REQUIRED_COLUMNS}")


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
        data: Input DataFrame.  Must contain the columns listed in REQUIRED_COLUMNS.
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
