"""Tests for eagar_tsai._core (integrand and compute_single_point)."""

from __future__ import annotations

import math

import pytest
import scipy.integrate

from eagar_tsai import BeamParameters, MaterialProperties, MeltPoolResult, SimulationDomain
from eagar_tsai._core import compute_single_point, eagar_tsai_integrand


class TestPythonIntegrand:
    """Tests for the pure-Python fallback integrand."""

    def test_positive_result(self) -> None:
        """Integrand is always positive for physical inputs."""
        assert eagar_tsai_integrand(1.0, 1.0, 0.0, 0.0, 0.5) > 0.0

    def test_symmetry_in_y(self) -> None:
        """Integrand is symmetric in y (even function)."""
        v1 = eagar_tsai_integrand(1.0, 0.5, 0.3, 0.0, 0.4)
        v2 = eagar_tsai_integrand(1.0, 0.5, -0.3, 0.0, 0.4)
        assert math.isclose(v1, v2, rel_tol=1e-12)

    def test_symmetry_in_z(self) -> None:
        """Integrand is symmetric in z."""
        v1 = eagar_tsai_integrand(1.0, 0.5, 0.0, 0.3, 0.4)
        v2 = eagar_tsai_integrand(1.0, 0.5, 0.0, -0.3, 0.4)
        assert math.isclose(v1, v2, rel_tol=1e-12)

    def test_integral_converges(self) -> None:
        """Integral from 0 to inf converges for typical parameters."""
        val, err = scipy.integrate.quad(eagar_tsai_integrand, 0.0, float("inf"), args=(0.5, 0.0, 0.0, 0.3))
        assert val > 0.0
        assert err < 1e-6

    def test_decays_far_from_source(self) -> None:
        """Integral value decreases for positions far from the heat source."""

        def integral_at(x: float) -> float:
            val, _ = scipy.integrate.quad(eagar_tsai_integrand, 0.0, float("inf"), args=(x, 0.0, 0.0, 0.3))
            return val

        assert integral_at(0.0) > integral_at(5.0) > integral_at(20.0)


@pytest.mark.integration
class TestCIntegrand:
    """Tests for the compiled C extension integrand."""

    def test_capsule_importable(self) -> None:
        """get_integrand_capsule() is importable and returns a non-None object."""
        from eagar_tsai._integrand_ext import get_integrand_capsule

        cap = get_integrand_capsule()
        assert cap is not None

    def test_lowlevelcallable_wraps(self) -> None:
        """LowLevelCallable wraps the capsule without error."""
        try:
            from scipy import LowLevelCallable
        except ImportError:
            from scipy.integrate import LowLevelCallable

        from eagar_tsai._integrand_ext import get_integrand_capsule

        llc = LowLevelCallable(get_integrand_capsule())
        assert llc is not None

    def test_result_matches_python_fallback(self) -> None:
        """C and Python integrands agree to within quad tolerance."""
        try:
            from scipy import LowLevelCallable
        except ImportError:
            from scipy.integrate import LowLevelCallable
        from scipy.integrate import quad

        from eagar_tsai._integrand_ext import get_integrand_capsule

        llc = LowLevelCallable(get_integrand_capsule())

        test_cases = [
            (0.5, 0.0, 0.0, 0.3),
            (1.0, 0.5, 0.2, 0.5),
            (2.0, 0.0, 1.0, 0.1),
        ]
        for x, y, z, p in test_cases:
            c_val, _ = quad(llc, 0.0, float("inf"), args=(x, y, z, p))
            py_val, _ = quad(eagar_tsai_integrand, 0.0, float("inf"), args=(x, y, z, p))
            assert math.isclose(c_val, py_val, rel_tol=1e-6), (
                f"Mismatch at (x={x}, y={y}, z={z}, p={p}): C={c_val:.6g}, Python={py_val:.6g}"
            )


class TestComputeSinglePointNonConvergence:
    """Tests for the domain-expansion failure path (no computation needed)."""

    def test_runtime_error_when_max_iters_zero(
        self,
        steel_beam: BeamParameters,
        steel_material: MaterialProperties,
        small_domain: SimulationDomain,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """RuntimeError is raised when _MAX_EXPANSION_ITERS is exhausted."""
        import eagar_tsai._core as core_module

        monkeypatch.setattr(core_module, "_MAX_EXPANSION_ITERS", 0)
        with pytest.raises(RuntimeError, match="Domain expansion did not converge"):
            compute_single_point(steel_beam, steel_material, small_domain)


@pytest.mark.slow
class TestComputeSinglePoint:
    """Integration tests for compute_single_point (marked slow)."""

    def test_returns_melt_pool_result(
        self,
        steel_beam: BeamParameters,
        steel_material: MaterialProperties,
        small_domain: SimulationDomain,
    ) -> None:
        """Return type is MeltPoolResult."""
        result = compute_single_point(steel_beam, steel_material, small_domain)
        assert isinstance(result, MeltPoolResult)

    def test_melt_pool_positive_dimensions(
        self,
        steel_beam: BeamParameters,
        steel_material: MaterialProperties,
        small_domain: SimulationDomain,
    ) -> None:
        """Dimensions are non-negative and peak temperature exceeds ambient."""
        result = compute_single_point(steel_beam, steel_material, small_domain)
        assert result.length >= 0.0
        assert result.width >= 0.0
        assert result.depth >= 0.0
        assert result.peak_temperature > 300.0

    def test_no_melt_pool_below_liquidus(
        self,
        small_domain: SimulationDomain,
    ) -> None:
        """Very low power produces no melt pool (all dimensions zero)."""
        beam = BeamParameters(
            beam_diameter=100e-6,
            power=1.0,  # extremely low
            velocity=5.0,
            absorptivity=0.35,
        )
        material = MaterialProperties(
            liquidus_temperature=5000.0,  # very high — impossible to melt
            thermal_conductivity=30.0,
            density=7800.0,
            specific_heat=700.0,
        )
        result = compute_single_point(beam, material, small_domain)
        assert result.length == 0.0
        assert result.width == 0.0
        assert result.depth == 0.0

    def test_peak_temperature_above_liquidus(
        self,
        steel_beam: BeamParameters,
        steel_material: MaterialProperties,
        small_domain: SimulationDomain,
    ) -> None:
        """If a melt pool exists, peak temperature must exceed liquidus."""
        result = compute_single_point(steel_beam, steel_material, small_domain)
        if result.length > 0.0:
            assert result.peak_temperature > steel_material.liquidus_temperature

    def test_default_domain_used_when_none(
        self,
        steel_beam: BeamParameters,
        steel_material: MaterialProperties,
    ) -> None:
        """Calling without a domain argument does not raise."""
        # This will use the default 1200x1200x1000 um, 1 um grid — very slow,
        # but we just check it doesn't error out. Marked slow for this reason.
        result = compute_single_point(steel_beam, steel_material)
        assert isinstance(result, MeltPoolResult)

    def test_micron_properties_consistent(
        self,
        steel_beam: BeamParameters,
        steel_material: MaterialProperties,
        small_domain: SimulationDomain,
    ) -> None:
        """Micrometre properties are exactly metres x 1e6."""
        result = compute_single_point(steel_beam, steel_material, small_domain)
        assert abs(result.length_um - result.length * 1e6) < 1e-12
        assert abs(result.width_um - result.width * 1e6) < 1e-12
        assert abs(result.depth_um - result.depth * 1e6) < 1e-12
