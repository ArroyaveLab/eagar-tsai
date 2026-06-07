"""Tests for eagar_tsai._api (DataFrame API and chunk processing)."""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from eagar_tsai import (
    SimulationDomain,
    TemperatureVolume,
    compute_melt_pool,
    compute_printability_map,
    compute_temperature_volumes,
)
from eagar_tsai._api import _REQUIRED_COLUMNS, _classify_defect, _process_chunk

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
    """Input validation tests (fast - no computation)."""

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
        """_REQUIRED_COLUMNS contains expected entries."""
        assert "velocity_m_s" in _REQUIRED_COLUMNS
        assert "power_w" in _REQUIRED_COLUMNS
        assert "absorptivity" in _REQUIRED_COLUMNS


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

        for col in _REQUIRED_COLUMNS:
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
        result = _process_chunk((0, _chunk(power_w=0.0), None, None, False))

        for col in _OUTPUT_COLUMNS:
            assert math.isnan(result.iloc[0][col]), f"Expected NaN for '{col}'"

    def test_error_row_preserves_input_columns(self) -> None:
        """Input column values survive even when the row errors."""
        result = _process_chunk((0, _chunk(power_w=0.0), None, None, False))

        assert "velocity_m_s" in result.columns
        assert result.iloc[0]["velocity_m_s"] == pytest.approx(0.5)

    def test_chunk_index_does_not_affect_nan_fill(self) -> None:
        """NaN fill is independent of the chunk index."""
        result_0 = _process_chunk((0, _chunk(power_w=0.0), None, None, False))
        result_7 = _process_chunk((7, _chunk(power_w=0.0), None, None, False))

        for col in _OUTPUT_COLUMNS:
            assert math.isnan(result_0.iloc[0][col])
            assert math.isnan(result_7.iloc[0][col])


@pytest.mark.slow
class TestProcessChunkCsvOutput:
    """output_dir triggers CSV writing to disk."""

    def test_csv_file_is_created(self, small_domain) -> None:
        """A CSV named ET_<chunk_idx>.csv is created inside output_dir."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            _process_chunk((3, _chunk(), small_domain, Path(tmp_dir), False))

            expected = Path(tmp_dir) / "ET_0003.csv"
            assert expected.exists(), f"Expected CSV not found: {expected}"

    def test_csv_content_matches_returned_dataframe(self, small_domain) -> None:
        """The CSV on disk contains the same columns and values as the return value."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            result = _process_chunk((0, _chunk(), small_domain, Path(tmp_dir), False))

            loaded = pd.read_csv(Path(tmp_dir) / "ET_0000.csv")

            assert list(loaded.columns) == list(result.columns)
            for col in _OUTPUT_COLUMNS:
                assert math.isclose(loaded.iloc[0][col], result.iloc[0][col], rel_tol=1e-9), f"Mismatch in column '{col}'"

    def test_output_dir_created_if_missing(self, small_domain) -> None:
        """mkdir(parents=True) means a nested path that does not yet exist is created."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            nested = Path(tmp_dir) / "nested" / "subdir"
            assert not nested.exists()

            _process_chunk((0, _chunk(), small_domain, nested, False))

            assert nested.is_dir()
            assert (nested / "ET_0000.csv").exists()

    def test_no_csv_when_output_dir_is_none(self, small_domain) -> None:
        """Passing output_dir=None leaves no files on disk."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            _process_chunk((0, _chunk(), small_domain, None, False))

            assert list(Path(tmp_dir).iterdir()) == [], "Unexpected files written"

    def test_mixed_chunk_only_valid_row_is_finite(self, small_domain) -> None:
        """In a chunk with one valid and one invalid row, only the invalid row gets NaN."""
        chunk = pd.DataFrame([_row(), _row(power_w=0.0)])

        with tempfile.TemporaryDirectory() as tmp_dir:
            result = _process_chunk((0, chunk, small_domain, Path(tmp_dir), False))

        for col in _OUTPUT_COLUMNS:
            assert math.isfinite(result.iloc[0][col]), f"Row 0 '{col}' should be finite"
            assert math.isnan(result.iloc[1][col]), f"Row 1 '{col}' should be NaN"


class TestClassifyDefect:
    """Unit tests for _classify_defect (no E-T computation - pure math)."""

    def test_no_melting_is_lof(self) -> None:
        """Zero depth or width returns lack of fusion."""
        defect, lof1, lof2, *_ = _classify_defect(0.0, 0.0, 0.0, 40.0, 90.0)
        assert defect == "lack_of_fusion"
        assert lof1 and lof2

    def test_lof1_shallow_depth(self) -> None:
        """Depth <= layer thickness triggers LOF1."""
        defect, lof1, lof2, ball1, ball2, kh1 = _classify_defect(200.0, 150.0, 20.0, 40.0, 90.0)
        assert lof1
        assert defect == "lack_of_fusion"

    def test_kh1_narrow_deep_pool(self) -> None:
        """Width/depth <= 2.5 triggers keyhole, taking priority over LOF."""
        defect, lof1, lof2, ball1, ball2, kh1 = _classify_defect(200.0, 50.0, 50.0, 40.0, 90.0)
        assert kh1
        assert defect == "keyhole"

    def test_ball1_elongated_pool(self) -> None:
        """Length/width >= 2.3 triggers balling when no other criterion fires."""
        defect, lof1, lof2, ball1, ball2, kh1 = _classify_defect(1500.0, 300.0, 100.0, 40.0, 90.0)
        assert ball1
        assert not kh1
        assert defect == "balling"

    def test_defect_free_window(self) -> None:
        """A well-proportioned melt pool with no criterion firing is defect-free."""
        defect, lof1, lof2, ball1, ball2, kh1 = _classify_defect(200.0, 250.0, 80.0, 40.0, 90.0)
        assert not kh1 and not lof1 and not lof2 and not ball1 and not ball2
        assert defect == "defect_free"


@pytest.mark.slow
class TestComputeTemperatureVolumes:
    """Integration tests for compute_temperature_volumes."""

    def test_returns_list_length_matches_input(self, volume_domain) -> None:
        """Returns one TemperatureVolume per input row."""
        df = pd.DataFrame([_row(), _row(power_w=250.0)])
        result = compute_temperature_volumes(df, domain=volume_domain, workers=1)
        assert len(result) == 2

    def test_items_are_temperature_volume(self, volume_domain) -> None:
        """Every returned item is a TemperatureVolume."""
        df = pd.DataFrame([_row()])
        result = compute_temperature_volumes(df, domain=volume_domain, workers=1)
        assert isinstance(result[0], TemperatureVolume)

    def test_raises_on_non_dataframe(self) -> None:
        """TypeError if input is not a DataFrame."""
        with pytest.raises(TypeError, match="pandas DataFrame"):
            compute_temperature_volumes({"not": "a dataframe"})  # ty: ignore[invalid-argument-type]

    def test_raises_on_missing_column(self) -> None:
        """ValueError lists the missing column."""
        df = pd.DataFrame([_row()])
        df = df.drop(columns=["power_w"])
        with pytest.raises(ValueError, match="power_w"):
            compute_temperature_volumes(df)


@pytest.mark.slow
class TestComputePrintabilityMap:
    """Integration tests for compute_printability_map."""

    def test_output_columns_present(self, steel_material, printability_domain, printability_params) -> None:
        """Result DataFrame contains all expected columns."""
        df = compute_printability_map(
            printability_params,
            steel_material,
            power_range=(100.0, 300.0),
            velocity_range=(0.5, 1.5),
            n_power=3,
            n_velocity=3,
            domain=printability_domain,
        )
        expected_cols = (
            "power_w",
            "velocity_m_s",
            "melt_length_um",
            "melt_width_um",
            "melt_depth_um",
            "defect",
            "lof1",
            "lof2",
            "ball1",
            "ball2",
            "kh1",
        )
        for col in expected_cols:
            assert col in df.columns, f"Missing column: {col}"

    def test_row_count_matches_grid(self, steel_material, printability_domain, printability_params) -> None:
        """Number of rows equals n_power * n_velocity."""
        df = compute_printability_map(
            printability_params,
            steel_material,
            power_range=(100.0, 300.0),
            velocity_range=(0.5, 1.5),
            n_power=3,
            n_velocity=4,
            domain=printability_domain,
        )
        assert len(df) == 12

    def test_defect_values_are_valid(self, steel_material, printability_domain, printability_params) -> None:
        """All defect labels are from the expected set."""
        df = compute_printability_map(
            printability_params,
            steel_material,
            power_range=(100.0, 300.0),
            velocity_range=(0.5, 1.5),
            n_power=3,
            n_velocity=3,
            domain=printability_domain,
        )
        assert set(df["defect"].unique()).issubset({"defect_free", "keyhole", "lack_of_fusion", "balling"})
