"""Core Eagar-Tsai solver: integrand and single-point melt pool computation.

The temperature field is computed by calling scipy.integrate.quad once per
grid point. When the compiled C extension is available, the integrand is
wrapped as a LowLevelCallable so QUADPACK calls it directly at C speed with
no Python overhead per evaluation. The pure-Python integrand is used as a
fallback when the extension cannot be imported.

Integrand (Sasha Rubenchik / LLNL 2015 reformulation of Eagar-Tsai 1983):

    f(t, x, y, z, p) = exp(-z^2/(4t) - (y^2 + (x-t)^2)/(4pt+1))
                       / ((4pt + 1) * sqrt(t))

where all coordinates are non-dimensionalised by the beam sigma parameter.
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

import numpy as np
from scipy.integrate import quad

from ._types import (
    _T0_K,
    BeamParameters,
    MaterialProperties,
    MeltPoolResult,
    SimulationDomain,
)

if TYPE_CHECKING:
    from collections.abc import Callable

__all__ = ["compute_single_point", "eagar_tsai_integrand"]

_logger = logging.getLogger(__name__)

_MAX_EXPANSION_ITERS: int = 20
"""Maximum number of domain-expansion iterations before giving up."""

_DEFAULT_DOMAIN = SimulationDomain()
"""Default simulation domain (1200 x 1200 x 1000 um, 1 um resolution)."""


def eagar_tsai_integrand(t: float, x: float, y: float, z: float, p: float) -> float:
    """Evaluate the Eagar-Tsai integrand at a single point.

    This function has the signature expected by scipy.integrate.quad
    when passed via the args keyword: quad(f, 0, inf, args=(x, y, z, p)).

    Args:
        t: Integration variable (dimensionless time), must be > 0.
        x: Non-dimensional x-coordinate (scan direction).
        y: Non-dimensional y-coordinate (cross-scan direction).
        z: Non-dimensional z-coordinate (depth, via sqrt(alpha * sigma / v)).
        p: Non-dimensional parameter alpha / (v * sigma).

    Returns:
        Integrand value at the given point.
    """
    denom_xy = 4.0 * p * t + 1.0
    z_exponent = -(z * z) / (4.0 * t)
    xy_exponent = -((y * y) + (x - t) ** 2) / denom_xy
    return math.exp(z_exponent + xy_exponent) / (denom_xy * math.sqrt(t))


def _get_integrand() -> Callable[..., float] | object:
    """Return the best available integrand for scipy.integrate.quad.

    Tries the compiled C extension first, returning a LowLevelCallable so
    QUADPACK can call the integrand at C speed with no Python overhead per
    evaluation. Falls back to the pure-Python implementation when the
    extension is unavailable.

    Returns:
        A LowLevelCallable (C extension) or the Python eagar_tsai_integrand.
    """
    try:
        try:
            from scipy import LowLevelCallable
        except ImportError:
            from scipy.integrate import LowLevelCallable  # ty: ignore[unresolved-import]

        from ._integrand_ext import get_integrand_capsule

        llc = LowLevelCallable(get_integrand_capsule())
        _logger.debug("Using compiled C integrand (LowLevelCallable).")
        return llc
    except ImportError:
        _logger.debug("C extension not available; using pure-Python fallback.")
        return eagar_tsai_integrand


_INTEGRAND = _get_integrand()
"""Cached integrand: LowLevelCallable (C extension) or Python fallback."""


def _build_grids(
    beam: BeamParameters,
    domain: SimulationDomain,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build 1-D coordinate arrays for the x-y and x-z temperature planes.

    Args:
        beam: Laser beam parameters.
        domain: Spatial domain specification.

    Returns:
        Tuple (xrange, yrange, zrange) of 1-D ndarrays in metres.
        xrange starts slightly behind the beam centre to capture the
        trailing melt pool; yrange runs from 0 to domain.y_length
        (half-domain, by symmetry); zrange runs from -domain.z_depth to 0.
    """
    delta = domain.spatial_resolution
    x_min = round(-1.5 * beam.beam_diameter, 9)  # slightly behind beam
    x_max = domain.x_length
    y_min = 0.0
    y_max = domain.y_length
    z_min = -domain.z_depth
    z_max = 0.0

    nx = int(np.round((x_max - x_min) / delta)) + 1
    ny = int(np.round((y_max - y_min) / delta)) + 1
    nz = int(np.round((z_max - z_min) / delta)) + 1

    return (
        np.linspace(x_min, x_max, nx),
        np.linspace(y_min, y_max, ny),
        np.linspace(z_min, z_max, nz),
    )


def _compute_temperature_planes(
    beam: BeamParameters,
    material: MaterialProperties,
    domain: SimulationDomain,
    integrand: Callable[..., float] | object,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the x-y (z=0) and x-z (y=0) temperature planes.

    Args:
        beam: Laser beam parameters.
        material: Material thermal properties.
        domain: Spatial domain specification.
        integrand: Callable accepted by scipy.integrate.quad — either the
            LowLevelCallable (C extension) or the Python fallback.

    Returns:
        Tuple (T_xy, T_xz) of 2-D ndarrays (temperature in Kelvin).
        T_xy has shape (ny, nx); T_xz has shape (nz, nx).
    """
    alpha = material.thermal_diffusivity
    sigma = beam.sigma
    velocity = beam.velocity
    absorptivity = beam.absorptivity
    power = beam.power
    conductivity = material.thermal_conductivity

    thermal_ratio = alpha / (velocity * sigma)
    temp_scale = (absorptivity * power) / (np.pi * (conductivity / alpha) * np.sqrt(np.pi * alpha * velocity * (sigma**3)))

    xrange, yrange, zrange = _build_grids(beam, domain)
    z_scale = np.sqrt((alpha * sigma) / velocity)

    T_xy = np.empty((yrange.size, xrange.size), dtype=np.float64)
    T_xz = np.empty((zrange.size, xrange.size), dtype=np.float64)

    for ix in range(xrange.size):
        x_nd = xrange[ix] / sigma
        for iy in range(yrange.size):
            val, _ = quad(integrand, 0.0, np.inf, args=(x_nd, yrange[iy] / sigma, 0.0, thermal_ratio))
            T_xy[iy, ix] = _T0_K + temp_scale * val
        for iz in range(zrange.size):
            val, _ = quad(integrand, 0.0, np.inf, args=(x_nd, 0.0, zrange[iz] / z_scale, thermal_ratio))
            T_xz[iz, ix] = _T0_K + temp_scale * val

    return T_xy, T_xz


def compute_single_point(
    beam: BeamParameters,
    material: MaterialProperties,
    domain: SimulationDomain | None = None,
) -> MeltPoolResult:
    """Compute melt pool dimensions for a single set of process parameters.

    The temperature field is evaluated on the x-y (z=0) and x-z (y=0)
    planes. Melt pool extents are extracted from the liquidus isotherm.
    If the melt pool touches a domain boundary the domain is expanded
    iteratively (up to _MAX_EXPANSION_ITERS times) and the computation is
    repeated on the enlarged grid.

    Args:
        beam: Laser beam and process parameters.
        material: Material thermal properties.
        domain: Spatial domain; defaults to _DEFAULT_DOMAIN
            (1200 x 1200 x 1000 um, 1 um).

    Returns:
        MeltPoolResult with melt pool geometry and temperature extremes.
        All lengths are zero if the peak temperature does not exceed the
        liquidus temperature.

    Raises:
        RuntimeError: If domain expansion does not converge within
            _MAX_EXPANSION_ITERS iterations.
    """
    if domain is None:
        domain = _DEFAULT_DOMAIN

    sigma_um = beam.sigma * 1.0e6  # metres -> um for domain expansion steps

    for iteration in range(_MAX_EXPANSION_ITERS):
        T_xy, T_xz = _compute_temperature_planes(beam, material, domain, _INTEGRAND)

        peak_T = float(np.amax(T_xy))
        min_T = float(np.amin(T_xy))

        if peak_T <= material.liquidus_temperature:
            return MeltPoolResult(
                length=0.0,
                width=0.0,
                depth=0.0,
                peak_temperature=peak_T,
                min_temperature=min_T,
            )

        xrange, yrange, zrange = _build_grids(beam, domain)
        t_melt = material.liquidus_temperature

        melt_x_mask = T_xy[0, :] > t_melt
        if not np.any(melt_x_mask):
            # Peak is above liquidus but the centreline row doesn't cross it —
            # shouldn't happen for a Gaussian source; treat as no melt pool.
            return MeltPoolResult(
                length=0.0,
                width=0.0,
                depth=0.0,
                peak_temperature=peak_T,
                min_temperature=min_T,
            )

        melt_x_indices = np.where(melt_x_mask)[0]
        melt_length = float(np.amax(xrange[melt_x_indices]) - np.amin(xrange[melt_x_indices]))

        y_half_length = 0.0
        z_length = 0.0
        for ix in melt_x_indices:
            melt_y_mask = T_xy[:, ix] > t_melt
            if np.any(melt_y_mask):
                melt_y_vals = yrange[melt_y_mask]
                tmp_y = float(np.amax(melt_y_vals) - np.amin(melt_y_vals))
                y_half_length = max(y_half_length, tmp_y)

            melt_z_mask = T_xz[:, ix] > t_melt
            if np.any(melt_z_mask):
                melt_z_vals = zrange[melt_z_mask]
                tmp_z = float(np.amax(melt_z_vals) - np.amin(melt_z_vals))
                z_length = max(z_length, tmp_z)

        x_max = domain.x_length
        y_max = domain.y_length
        z_span = domain.z_depth  # positive

        needs_expansion_x = np.isclose(float(np.amax(xrange[melt_x_indices])), x_max)
        needs_expansion_y = np.isclose(y_half_length, y_max)
        needs_expansion_z = np.isclose(z_length, z_span)

        if needs_expansion_x:
            _logger.info(
                "Iteration %d: x domain too small (%.1f um), expanding by %.1f um.",
                iteration,
                domain.x_length_um,
                sigma_um,
            )
            domain = domain.expanded(dx_um=sigma_um)
            continue

        if needs_expansion_y:
            _logger.info(
                "Iteration %d: y domain too small (%.1f um), expanding by %.1f um.",
                iteration,
                domain.y_length_um,
                sigma_um,
            )
            domain = domain.expanded(dy_um=sigma_um)
            continue

        if needs_expansion_z:
            _logger.info(
                "Iteration %d: z domain too small (%.1f um), expanding by %.1f um.",
                iteration,
                domain.z_depth_um,
                sigma_um,
            )
            domain = domain.expanded(dz_um=sigma_um)
            continue

        return MeltPoolResult(
            length=melt_length,
            width=y_half_length * 2.0,  # half-domain symmetry
            depth=z_length,
            peak_temperature=peak_T,
            min_temperature=min_T,
        )

    raise RuntimeError(
        f"Domain expansion did not converge after {_MAX_EXPANSION_ITERS} iterations. Consider increasing the default domain size."
    )
