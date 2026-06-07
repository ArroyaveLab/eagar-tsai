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

import concurrent.futures
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
    TemperatureField,
    TemperatureVolume,
)

if TYPE_CHECKING:
    from collections.abc import Callable

__all__ = ["compute_single_point", "compute_temperature_volume", "eagar_tsai_integrand"]

_logger = logging.getLogger(__name__)

_MAX_EXPANSION_ITERS: int = 20
"""Maximum number of domain-expansion iterations before giving up."""

"""Default simulation domain (1200 x 1200 x 1000 um, 1 um resolution)."""


def eagar_tsai_integrand(t: float, x: float, y: float, z: float, p: float) -> float:
    """Evaluate the Eagar-Tsai integrand at a single point.

    This function has the signature expected by scipy.integrate.quad
    when passed via the args keyword: quad(f, 0, inf, args=(x, y, z, p)).

    Args:
        t: Integration variable (dimensionless time), must be > 0.
        x: Non-dimensional x-coordinate, positive in the trailing wake direction
            (opposite to the scan direction).
        y: Non-dimensional y-coordinate (cross-scan direction).
        z: Non-dimensional z-coordinate (depth, scaled by sqrt(alpha * sigma / v)).
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
        xrange starts slightly ahead of the beam center in the scan direction
        (x < 0); the trailing melt pool extends into positive x.
        yrange runs from 0 to domain.y_length
        (half-domain, by symmetry); zrange runs from -domain.z_depth to 0.
    """
    delta = domain.spatial_resolution
    x_min = round(-1.5 * beam.beam_diameter, 9)  # small margin ahead of beam in scan direction
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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the x-y (z=0) and x-z (y=0) temperature planes.

    Args:
        beam: Laser beam parameters.
        material: Material thermal properties.
        domain: Spatial domain specification.
        integrand: Callable accepted by scipy.integrate.quad — either the
            LowLevelCallable (C extension) or the Python fallback.

    Returns:
        Tuple (T_xy, T_xz, xrange, yrange, zrange).
        T_xy has shape (ny, nx); T_xz has shape (nz, nx).
        Coordinate arrays are in metres.
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

    return T_xy, T_xz, xrange, yrange, zrange


def compute_single_point(
    beam: BeamParameters,
    material: MaterialProperties,
    domain: SimulationDomain | None = None,
    *,
    full_field: bool = True,
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
        domain: Spatial domain; defaults to 1200 x 1200 x 1000 um, 1 um.
        full_field: When ``True`` (default), T_xy and T_xz are computed for
            every grid point and stored in the returned ``TemperatureField``.
            When ``False``, T_xz is computed only at the x-indices that fall
            inside the melt pool, which is significantly faster when the melt
            pool spans a small fraction of the domain. The TemperatureField
            is still returned but T_xz is filled with the ambient temperature
            outside the melt x-range, so it should not be used for plotting.

    Returns:
        A ``MeltPoolResult`` containing melt pool dimensions, temperature
        extremes, and the full ``TemperatureField``. All lengths are zero
        if the peak temperature does not exceed the liquidus temperature.

    Raises:
        RuntimeError: If domain expansion does not converge within _MAX_EXPANSION_ITERS iterations.
    """
    if domain is None:
        domain = SimulationDomain()

    sigma_um = beam.sigma * 1.0e6  # metres -> um for domain expansion steps

    # Thermal scaling factors depend only on beam/material, not on the domain,
    # so they are computed once outside the expansion loop.
    alpha = material.thermal_diffusivity
    sigma = beam.sigma
    velocity = beam.velocity
    thermal_ratio = alpha / (velocity * sigma)
    temp_scale = (beam.absorptivity * beam.power) / (
        np.pi * (material.thermal_conductivity / alpha) * np.sqrt(np.pi * alpha * velocity * (sigma**3))
    )
    z_scale = np.sqrt((alpha * sigma) / velocity)

    for iteration in range(_MAX_EXPANSION_ITERS):
        if full_field:
            T_xy, T_xz, xrange, yrange, zrange = _compute_temperature_planes(beam, material, domain, _INTEGRAND)
        else:
            xrange, yrange, zrange = _build_grids(beam, domain)
            T_xy = np.empty((yrange.size, xrange.size), dtype=np.float64)
            for ix in range(xrange.size):
                x_nd = xrange[ix] / sigma
                for iy in range(yrange.size):
                    val, _ = quad(_INTEGRAND, 0.0, np.inf, args=(x_nd, yrange[iy] / sigma, 0.0, thermal_ratio))
                    T_xy[iy, ix] = _T0_K + temp_scale * val

        peak_T = float(np.amax(T_xy))
        min_T = float(np.amin(T_xy))

        if peak_T <= material.liquidus_temperature:
            if not full_field:
                T_xz = np.full((zrange.size, xrange.size), _T0_K, dtype=np.float64)
            tf = TemperatureField(
                T_xy=T_xy,
                T_xz=T_xz,
                x_range_m=xrange,
                y_range_m=yrange,
                z_range_m=zrange,
                liquidus_temperature_k=material.liquidus_temperature,
                melt_width_m=0.0,
                melt_depth_m=0.0,
            )
            return MeltPoolResult(
                length=0.0,
                width=0.0,
                depth=0.0,
                peak_temperature=peak_T,
                min_temperature=min_T,
                temperature_field=tf,
            )

        t_melt = material.liquidus_temperature

        melt_x_mask = T_xy[0, :] > t_melt
        if not np.any(melt_x_mask):
            # Peak is above liquidus but the centerline row doesn't cross it
            # shouldn't happen for a Gaussian source; treat as no melt pool.
            if not full_field:
                T_xz = np.full((zrange.size, xrange.size), _T0_K, dtype=np.float64)
            tf = TemperatureField(
                T_xy=T_xy,
                T_xz=T_xz,
                x_range_m=xrange,
                y_range_m=yrange,
                z_range_m=zrange,
                liquidus_temperature_k=material.liquidus_temperature,
                melt_width_m=0.0,
                melt_depth_m=0.0,
            )
            return MeltPoolResult(
                length=0.0,
                width=0.0,
                depth=0.0,
                peak_temperature=peak_T,
                min_temperature=min_T,
                temperature_field=tf,
            )

        melt_x_indices = np.where(melt_x_mask)[0]
        melt_length = float(np.amax(xrange[melt_x_indices]) - np.amin(xrange[melt_x_indices]))

        x_max = domain.x_length
        y_max = domain.y_length
        z_span = domain.z_depth

        if full_field:
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

        else:
            # Fast path: check x expansion before paying the cost of T_xz.
            needs_expansion_x = np.isclose(float(np.amax(xrange[melt_x_indices])), x_max)
            if needs_expansion_x:
                _logger.info(
                    "Iteration %d: x domain too small (%.1f um), expanding by %.1f um.",
                    iteration,
                    domain.x_length_um,
                    sigma_um,
                )
                domain = domain.expanded(dx_um=sigma_um)
                continue

            # Extract y extent from T_xy (already computed, no new quad calls).
            y_half_length = 0.0
            for ix in melt_x_indices:
                melt_y_mask = T_xy[:, ix] > t_melt
                if np.any(melt_y_mask):
                    melt_y_vals = yrange[melt_y_mask]
                    y_half_length = max(y_half_length, float(np.amax(melt_y_vals) - np.amin(melt_y_vals)))

            needs_expansion_y = np.isclose(y_half_length, y_max)
            if needs_expansion_y:
                _logger.info(
                    "Iteration %d: y domain too small (%.1f um), expanding by %.1f um.",
                    iteration,
                    domain.y_length_um,
                    sigma_um,
                )
                domain = domain.expanded(dy_um=sigma_um)
                continue

            # Compute T_xz only for x-indices inside the melt pool.
            T_xz = np.full((zrange.size, xrange.size), _T0_K, dtype=np.float64)
            for ix in melt_x_indices:
                x_nd = xrange[ix] / sigma
                for iz in range(zrange.size):
                    val, _ = quad(_INTEGRAND, 0.0, np.inf, args=(x_nd, 0.0, zrange[iz] / z_scale, thermal_ratio))
                    T_xz[iz, ix] = _T0_K + temp_scale * val

            z_length = 0.0
            for ix in melt_x_indices:
                melt_z_mask = T_xz[:, ix] > t_melt
                if np.any(melt_z_mask):
                    melt_z_vals = zrange[melt_z_mask]
                    z_length = max(z_length, float(np.amax(melt_z_vals) - np.amin(melt_z_vals)))

            needs_expansion_z = np.isclose(z_length, z_span)
            if needs_expansion_z:
                _logger.info(
                    "Iteration %d: z domain too small (%.1f um), expanding by %.1f um.",
                    iteration,
                    domain.z_depth_um,
                    sigma_um,
                )
                domain = domain.expanded(dz_um=sigma_um)
                continue

        melt_width = y_half_length * 2.0  # half-domain symmetry
        tf = TemperatureField(
            T_xy=T_xy,
            T_xz=T_xz,
            x_range_m=xrange,
            y_range_m=yrange,
            z_range_m=zrange,
            liquidus_temperature_k=material.liquidus_temperature,
            melt_width_m=melt_width,
            melt_depth_m=z_length,
        )
        return MeltPoolResult(
            length=melt_length,
            width=melt_width,
            depth=z_length,
            peak_temperature=peak_T,
            min_temperature=min_T,
            temperature_field=tf,
        )

    raise RuntimeError(
        f"Domain expansion did not converge after {_MAX_EXPANSION_ITERS} iterations. Consider increasing the default domain size."
    )


def _compute_volume_x_slice(
    args: tuple[int, int, np.ndarray, np.ndarray, np.ndarray, float, float, float, float],
) -> np.ndarray:
    """Compute temperature for a contiguous range of x-slices in the 3-D volume.

    Top-level module function so it can be pickled for ``ProcessPoolExecutor``.

    Args:
        args: Tuple of ``(ix_start, ix_end, x_range, y_range, z_range,
            sigma, z_scale, thermal_ratio, temp_scale)``.
            Coordinate arrays must be in metres.

    Returns:
        Array of shape ``(ix_end - ix_start, ny, nz)`` in Kelvin.
    """
    ix_start, ix_end, x_range, y_range, z_range, sigma, z_scale, thermal_ratio, temp_scale = args
    ny = y_range.size
    nz = z_range.size
    T_slice = np.empty((ix_end - ix_start, ny, nz), dtype=np.float64)
    for i, ix in enumerate(range(ix_start, ix_end)):
        x_nd = x_range[ix] / sigma
        for iy in range(ny):
            y_nd = y_range[iy] / sigma
            for iz in range(nz):
                val, _ = quad(_INTEGRAND, 0.0, np.inf, args=(x_nd, y_nd, z_range[iz] / z_scale, thermal_ratio))
                T_slice[i, iy, iz] = _T0_K + temp_scale * val
    return T_slice


def compute_temperature_volume(
    beam: BeamParameters,
    material: MaterialProperties,
    domain: SimulationDomain | None = None,
    *,
    workers: int | None = None,
    chunk_size: int = 10,
) -> TemperatureVolume:
    """Compute the full 3-D temperature volume for a single set of process parameters.

    The spatial domain is auto-sized by first running the standard 2-D melt
    pool computation (``compute_single_point``), which iteratively expands the
    domain until the melt pool fits within all boundaries. The 3-D volume is
    then evaluated on that converged domain.

    Computation is parallelised over contiguous x-index slices using
    ``concurrent.futures.ProcessPoolExecutor`` when ``workers > 1``.

    Args:
        beam: Laser beam and process parameters.
        material: Material thermal properties.
        domain: Starting simulation domain for the auto-sizing step.
            When ``None``, the default ``SimulationDomain()`` is used.
        workers: Worker processes for parallel x-slice computation.
            ``None`` or ``1`` runs serially. ``-1`` uses all available cores.
        chunk_size: Number of x-index slices dispatched to each worker task.
            Larger values reduce multiprocessing overhead; smaller values
            improve load balancing when slice cost varies.

    Returns:
        A ``TemperatureVolume`` containing the 3-D temperature array of shape
        ``(nx, ny, nz)``, coordinate arrays in metres, the liquidus temperature,
        and the ``MeltPoolResult`` from the auto-sizing step.

    Raises:
        RuntimeError: If domain expansion does not converge within
            ``_MAX_EXPANSION_ITERS`` iterations.

    Examples:
        ```python
        from eagar_tsai import BeamParameters, MaterialProperties, compute_temperature_volume

        beam = BeamParameters(beam_diameter=80e-6, power=250.0, velocity=0.5, absorptivity=0.59)
        mat = MaterialProperties(liquidus_temperature=3455.0, thermal_conductivity=23.75,
                                 density=18038.9, specific_heat=251.6)
        vol = compute_temperature_volume(beam, mat, workers=-1)
        vol.export_vti("temperature_volume.vti")
        vol.plot_3d()
        ```
    """
    result = compute_single_point(beam, material, domain, full_field=False)

    tf = result.temperature_field
    xrange = tf.x_range_m
    yrange = tf.y_range_m
    zrange = tf.z_range_m

    alpha = material.thermal_diffusivity
    sigma = beam.sigma
    velocity = beam.velocity
    thermal_ratio = alpha / (velocity * sigma)
    temp_scale = (beam.absorptivity * beam.power) / (
        np.pi * (material.thermal_conductivity / alpha) * np.sqrt(np.pi * alpha * velocity * (sigma**3))
    )
    z_scale = np.sqrt((alpha * sigma) / velocity)

    nx = xrange.size
    x_starts = list(range(0, nx, chunk_size))
    x_ends = [min(s + chunk_size, nx) for s in x_starts]
    slice_args = [
        (s, e, xrange, yrange, zrange, sigma, z_scale, thermal_ratio, temp_scale) for s, e in zip(x_starts, x_ends, strict=True)
    ]

    if workers is None or workers == 1:
        slices = [_compute_volume_x_slice(a) for a in slice_args]
    else:
        max_workers = None if workers == -1 else workers
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            slices = list(executor.map(_compute_volume_x_slice, slice_args))

    T_xyz = np.concatenate(slices, axis=0)  # (nx, ny, nz)

    return TemperatureVolume(
        T_xyz=T_xyz,
        x_range_m=xrange,
        y_range_m=yrange,
        z_range_m=zrange,
        liquidus_temperature_k=material.liquidus_temperature,
        result=result,
    )
