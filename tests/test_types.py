"""Tests for frozen dataclasses in eagar_tsai._types."""

from __future__ import annotations

import math

import pytest

from eagar_tsai import BeamParameters, MaterialProperties, MeltPoolResult, SimulationDomain
from eagar_tsai._types import _T0_K


class TestBeamParameters:
    """Tests for BeamParameters."""

    def test_sigma_derived_correctly(self) -> None:
        """Sigma = sqrt(2) * beam_diameter / 2."""
        beam = BeamParameters(
            beam_diameter=100e-6,
            power=200.0,
            velocity=0.5,
            absorptivity=0.35,
        )
        expected = math.sqrt(2.0) * 50e-6
        assert math.isclose(beam.sigma, expected, rel_tol=1e-12)

    def test_frozen(self) -> None:
        """Dataclass fields cannot be mutated."""
        beam = BeamParameters(beam_diameter=100e-6, power=200.0, velocity=0.5, absorptivity=0.35)
        with pytest.raises((AttributeError, TypeError)):
            beam.power = 300.0  # type: ignore[misc]

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"beam_diameter": -1e-6, "power": 200.0, "velocity": 0.5, "absorptivity": 0.35},
            {"beam_diameter": 100e-6, "power": 0.0, "velocity": 0.5, "absorptivity": 0.35},
            {"beam_diameter": 100e-6, "power": 200.0, "velocity": -1.0, "absorptivity": 0.35},
            {"beam_diameter": 100e-6, "power": 200.0, "velocity": 0.5, "absorptivity": 0.0},
            {"beam_diameter": 100e-6, "power": 200.0, "velocity": 0.5, "absorptivity": 1.5},
        ],
    )
    def test_invalid_values(self, kwargs: dict) -> None:
        """Invalid physical values raise ValueError."""
        with pytest.raises(ValueError):
            BeamParameters(**kwargs)


class TestMaterialProperties:
    """Tests for MaterialProperties."""

    def test_thermal_diffusivity(self) -> None:
        """thermal_diffusivity = k / (rho * cp)."""
        mat = MaterialProperties(
            liquidus_temperature=1700.0,
            thermal_conductivity=30.0,
            density=7800.0,
            specific_heat=700.0,
        )
        expected = 30.0 / (7800.0 * 700.0)
        assert math.isclose(mat.thermal_diffusivity, expected, rel_tol=1e-12)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {
                "liquidus_temperature": 0.0,
                "thermal_conductivity": 30.0,
                "density": 7800.0,
                "specific_heat": 700.0,
            },
            {
                "liquidus_temperature": 1700.0,
                "thermal_conductivity": -1.0,
                "density": 7800.0,
                "specific_heat": 700.0,
            },
            {
                "liquidus_temperature": 1700.0,
                "thermal_conductivity": 30.0,
                "density": 0.0,
                "specific_heat": 700.0,
            },
            {
                "liquidus_temperature": 1700.0,
                "thermal_conductivity": 30.0,
                "density": 7800.0,
                "specific_heat": -1.0,
            },
        ],
    )
    def test_invalid_values(self, kwargs: dict) -> None:
        """Invalid physical values raise ValueError."""
        with pytest.raises(ValueError):
            MaterialProperties(**kwargs)


class TestSimulationDomain:
    """Tests for SimulationDomain."""

    def test_metre_properties(self) -> None:
        """Property accessors return correct SI values."""
        d = SimulationDomain(x_length_um=1200.0, y_length_um=800.0, z_depth_um=500.0, spatial_resolution_um=2.0)
        assert math.isclose(d.x_length, 1200e-6, rel_tol=1e-12)
        assert math.isclose(d.y_length, 800e-6, rel_tol=1e-12)
        assert math.isclose(d.z_depth, 500e-6, rel_tol=1e-12)
        assert math.isclose(d.spatial_resolution, 2e-6, rel_tol=1e-12)

    def test_expanded_returns_new_instance(self) -> None:
        """expanded() returns a new object without mutating the original."""
        d = SimulationDomain(x_length_um=1000.0, y_length_um=1000.0, z_depth_um=800.0)
        d2 = d.expanded(dx_um=100.0, dy_um=50.0, dz_um=25.0)
        assert d2 is not d
        assert math.isclose(d2.x_length_um, 1100.0)
        assert math.isclose(d2.y_length_um, 1050.0)
        assert math.isclose(d2.z_depth_um, 825.0)
        assert math.isclose(d.x_length_um, 1000.0)

    def test_frozen(self) -> None:
        """SimulationDomain is immutable."""
        d = SimulationDomain()
        with pytest.raises((AttributeError, TypeError)):
            d.x_length_um = 500.0  # type: ignore[misc]

    def test_invalid_negative_dimension(self) -> None:
        """Non-positive dimension raises ValueError."""
        with pytest.raises(ValueError):
            SimulationDomain(x_length_um=-100.0)


class TestMeltPoolResult:
    """Tests for MeltPoolResult."""

    def test_micron_properties(self) -> None:
        """Micrometre properties are metres x 1e6."""
        r = MeltPoolResult(
            length=500e-6,
            width=300e-6,
            depth=150e-6,
            peak_temperature=2000.0,
            min_temperature=300.0,
        )
        assert math.isclose(r.length_um, 500.0)
        assert math.isclose(r.width_um, 300.0)
        assert math.isclose(r.depth_um, 150.0)


class TestConstants:
    """Tests for module-level constants."""

    def test_ambient_temperature(self) -> None:
        """Ambient temperature constant is 300 K."""
        assert _T0_K == 300.0
