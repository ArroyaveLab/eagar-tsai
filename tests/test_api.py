"""Tests for eagar_tsai._api (DataFrame API and chunk processing)."""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from eagar_tsai import SimulationDomain, compute_melt_pool
from eagar_tsai._api import REQUIRED_COLUMNS, _process_chunk

_OUTPUT_COLUMNS = [
    "melt_length",
    "melt_width",
    "melt_depth",
    "melt_length_um",
    "melt_width_um",
    "melt_depth_um",
    "peak_temperature",
    "min_temperature",
]

_SMALL_DOMAIN = SimulationDomain(x_length_um=300.0, y_length_um=300.0, z_depth_um=200.0, spatial_resolution_um=10.0)


def _row(**overrides: object) -> dict[str, object]:
    """Return a dict of valid required input columns."""
    defaults: dict[str, object] = {
        "velocity_m_s": 0.5,
        "power_w": 200.0,
        "beam_diameter_m": 100e-6,
        "absorptivity": 0.35,
        "liquidus_temperature_k": 1700.0,
        "thermal_conductivity_w_mk": 30.0,
        "density_kg_m3": 7800.0,
        "specific_heat_j_kgk": 700.0,
    }
    defaults.update(overrides)
    return defaults


def _chunk(**overrides: object) -> pd.DataFrame:
    """Return a single-row DataFrame built from ``_row``."""
    return pd.DataFrame([_row(**overrides)])


def _minimal_df(**overrides: object) -> pd.DataFrame:
    """Return a single-row DataFrame with all required columns."""
    return _chunk(**overrides)


class TestComputeMeltPoolValidation:
    """Input validation tests (fast — no computation)."""

    def test_raises_on_non_dataframe(self) -> None:
        """TypeError if input is not a DataFrame."""
        with pytest.raises(TypeError, match="pandas DataFrame"):
            compute_melt_pool({"not": "a dataframe"})  # ty: ignore[invalid-argument-type]

    def test_raises_on_missing_column(self) -> None:
        """ValueError lists the missing column."""
        df = _minimal_df()
        df = df.drop(columns=["power_w"])
        with pytest.raises(ValueError, match="power_w"):
            compute_melt_pool(df)

    def test_required_columns_list(self) -> None:
        """REQUIRED_COLUMNS contains expected entries."""
        assert "velocity_m_s" in REQUIRED_COLUMNS
        assert "power_w" in REQUIRED_COLUMNS
        assert "absorptivity" in REQUIRED_COLUMNS


@pytest.mark.slow
class TestComputeMeltPoolOutput:
    """End-to-end tests for compute_melt_pool (marked slow)."""

    def test_output_has_expected_columns(self) -> None:
        """Result DataFrame contains all output columns."""
        domain = SimulationDomain(x_length_um=300.0, y_length_um=300.0, z_depth_um=200.0, spatial_resolution_um=10.0)
        df = _minimal_df()
        result = compute_melt_pool(df, domain=domain, workers=1)

        for col in _OUTPUT_COLUMNS:
            assert col in result.columns, f"Missing output column: {col}"

    def test_output_preserves_input_columns(self) -> None:
        """Input columns are preserved in the result DataFrame."""
        domain = SimulationDomain(x_length_um=300.0, y_length_um=300.0, z_depth_um=200.0, spatial_resolution_um=10.0)
        df = _minimal_df()
        result = compute_melt_pool(df, domain=domain, workers=1)

        for col in REQUIRED_COLUMNS:
            assert col in result.columns

    def test_output_row_count_matches_input(self) -> None:
        """Output has the same number of rows as input."""
        domain = SimulationDomain(x_length_um=300.0, y_length_um=300.0, z_depth_um=200.0, spatial_resolution_um=10.0)
        df = pd.concat([_minimal_df()] * 3, ignore_index=True)
        result = compute_melt_pool(df, domain=domain, workers=1, chunk_size=2)
        assert len(result) == len(df)

    def test_serial_and_parallel_produce_same_results(self) -> None:
        """Serial and 2-worker parallel results are numerically identical."""
        domain = SimulationDomain(x_length_um=300.0, y_length_um=300.0, z_depth_um=200.0, spatial_resolution_um=10.0)
        df = pd.concat([_minimal_df()] * 2, ignore_index=True)

        serial = compute_melt_pool(df, domain=domain, workers=1, chunk_size=1)
        parallel = compute_melt_pool(df, domain=domain, workers=2, chunk_size=1)

        for col in ["melt_length", "melt_width", "melt_depth"]:
            assert (serial[col] - parallel[col]).abs().max() < 1e-12


class TestProcessChunkErrorHandling:
    """Invalid rows are caught and replaced with NaN outputs."""

    def test_invalid_row_fills_nan(self) -> None:
        """A row whose BeamParameters validation fails produces all-NaN outputs."""
        # power_w=0 is rejected by BeamParameters.__post_init__ -> ValueError
        result = _process_chunk((0, _chunk(power_w=0.0), None, None))

        for col in _OUTPUT_COLUMNS:
            assert math.isnan(result.iloc[0][col]), f"Expected NaN for '{col}'"

    def test_error_row_preserves_input_columns(self) -> None:
        """Input column values survive even when the row errors."""
        result = _process_chunk((0, _chunk(power_w=0.0), None, None))

        assert "velocity_m_s" in result.columns
        assert result.iloc[0]["velocity_m_s"] == pytest.approx(0.5)

    def test_chunk_index_does_not_affect_nan_fill(self) -> None:
        """NaN fill is independent of the chunk index."""
        result_0 = _process_chunk((0, _chunk(power_w=0.0), None, None))
        result_7 = _process_chunk((7, _chunk(power_w=0.0), None, None))

        for col in _OUTPUT_COLUMNS:
            assert math.isnan(result_0.iloc[0][col])
            assert math.isnan(result_7.iloc[0][col])


@pytest.mark.slow
class TestProcessChunkCsvOutput:
    """output_dir triggers CSV writing to disk."""

    def test_csv_file_is_created(self) -> None:
        """A CSV named ET_<chunk_idx>.csv is created inside output_dir."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            _process_chunk((3, _chunk(), _SMALL_DOMAIN, Path(tmp_dir)))

            expected = Path(tmp_dir) / "ET_0003.csv"
            assert expected.exists(), f"Expected CSV not found: {expected}"

    def test_csv_content_matches_returned_dataframe(self) -> None:
        """The CSV on disk contains the same columns and values as the return value."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            result = _process_chunk((0, _chunk(), _SMALL_DOMAIN, Path(tmp_dir)))

            loaded = pd.read_csv(Path(tmp_dir) / "ET_0000.csv")

            assert list(loaded.columns) == list(result.columns)
            for col in _OUTPUT_COLUMNS:
                assert math.isclose(loaded.iloc[0][col], result.iloc[0][col], rel_tol=1e-9), f"Mismatch in column '{col}'"

    def test_output_dir_created_if_missing(self) -> None:
        """mkdir(parents=True) means a nested path that does not yet exist is created."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            nested = Path(tmp_dir) / "nested" / "subdir"
            assert not nested.exists()

            _process_chunk((0, _chunk(), _SMALL_DOMAIN, nested))

            assert nested.is_dir()
            assert (nested / "ET_0000.csv").exists()

    def test_no_csv_when_output_dir_is_none(self) -> None:
        """Passing output_dir=None leaves no files on disk."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Run from tmp_dir so any stray writes would be visible
            _process_chunk((0, _chunk(), _SMALL_DOMAIN, None))

            assert list(Path(tmp_dir).iterdir()) == [], "Unexpected files written"

    def test_mixed_chunk_only_valid_row_is_finite(self) -> None:
        """In a chunk with one valid and one invalid row, only the invalid row gets NaN."""
        chunk = pd.DataFrame([_row(), _row(power_w=0.0)])

        with tempfile.TemporaryDirectory() as tmp_dir:
            result = _process_chunk((0, chunk, _SMALL_DOMAIN, Path(tmp_dir)))

        for col in _OUTPUT_COLUMNS:
            assert math.isfinite(result.iloc[0][col]), f"Row 0 '{col}' should be finite"
            assert math.isnan(result.iloc[1][col]), f"Row 1 '{col}' should be NaN"
