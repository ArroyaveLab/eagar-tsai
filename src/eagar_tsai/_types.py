"""Frozen dataclasses for Eagar–Tsai model parameters and results."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

__all__ = [
    "BeamParameters",
    "MaterialProperties",
    "SimulationDomain",
    "MeltPoolResult",
]

_T0_K: float = 300.0
"""Ambient temperature in Kelvin (hard-coded per the original model)."""

_UM_TO_M: float = 1.0e-6
"""Conversion factor from micrometres to metres."""


@dataclass(frozen=True)
class BeamParameters:
    """Laser beam and process parameters.

    Attributes:
        beam_diameter: Beam diameter (2*sigma) in metres.
        power: Laser power in Watts.
        velocity: Scan velocity in m/s.
        absorptivity: Absorptivity (dimensionless, must be in (0, 1]).
        sigma: Derived Gaussian width parameter sqrt(2) * (beam_diameter / 2).
    """

    beam_diameter: float
    power: float
    velocity: float
    absorptivity: float
    sigma: float = field(init=False)

    def __post_init__(self) -> None:
        """Validate physical constraints and compute derived sigma."""
        if self.beam_diameter <= 0.0:
            raise ValueError(f"beam_diameter must be positive, got {self.beam_diameter}")
        if self.power <= 0.0:
            raise ValueError(f"power must be positive, got {self.power}")
        if self.velocity <= 0.0:
            raise ValueError(f"velocity must be positive, got {self.velocity}")
        if not (0.0 < self.absorptivity <= 1.0):
            raise ValueError(f"absorptivity must be in (0, 1], got {self.absorptivity}")
        # sigma = sqrt(2) * (beam_diameter / 2) — 1/e^2 Gaussian half-width
        object.__setattr__(self, "sigma", math.sqrt(2.0) * (self.beam_diameter / 2.0))


@dataclass(frozen=True)
class MaterialProperties:
    """Material thermal properties evaluated at the liquidus temperature.

    Attributes:
        liquidus_temperature: Liquidus temperature in Kelvin.
        thermal_conductivity: Thermal conductivity in W/(m K).
        density: Density in kg/m^3.
        specific_heat: Specific heat capacity in J/(kg K).
    """

    liquidus_temperature: float
    thermal_conductivity: float
    density: float
    specific_heat: float

    def __post_init__(self) -> None:
        """Validate physical constraints."""
        if self.liquidus_temperature <= 0.0:
            raise ValueError(f"liquidus_temperature must be positive, got {self.liquidus_temperature}")
        if self.thermal_conductivity <= 0.0:
            raise ValueError(f"thermal_conductivity must be positive, got {self.thermal_conductivity}")
        if self.density <= 0.0:
            raise ValueError(f"density must be positive, got {self.density}")
        if self.specific_heat <= 0.0:
            raise ValueError(f"specific_heat must be positive, got {self.specific_heat}")

    @property
    def thermal_diffusivity(self) -> float:
        """Thermal diffusivity alpha = k / (rho * cp) in m^2/s."""
        return self.thermal_conductivity / (self.density * self.specific_heat)


@dataclass(frozen=True)
class SimulationDomain:
    """Spatial domain for numerical temperature field evaluation.

    All _um attributes are in micrometres; corresponding properties
    (without the suffix) return the equivalent value in metres.

    Attributes:
        x_length_um: Domain length along x (scan direction) in um.
        y_length_um: Domain half-width along y in um.
        z_depth_um: Domain depth along z in um.
        spatial_resolution_um: Grid spacing in um.
    """

    x_length_um: float = 1200.0
    y_length_um: float = 1200.0
    z_depth_um: float = 1000.0
    spatial_resolution_um: float = 1.0

    def __post_init__(self) -> None:
        """Validate that all domain dimensions and resolution are positive."""
        for attr in ("x_length_um", "y_length_um", "z_depth_um", "spatial_resolution_um"):
            val = getattr(self, attr)
            if val <= 0.0:
                raise ValueError(f"{attr} must be positive, got {val}")

    @property
    def x_length(self) -> float:
        """Domain length along x in metres."""
        return self.x_length_um * _UM_TO_M

    @property
    def y_length(self) -> float:
        """Domain half-width along y in metres."""
        return self.y_length_um * _UM_TO_M

    @property
    def z_depth(self) -> float:
        """Domain depth along z in metres."""
        return self.z_depth_um * _UM_TO_M

    @property
    def spatial_resolution(self) -> float:
        """Grid spacing in metres."""
        return self.spatial_resolution_um * _UM_TO_M

    def expanded(
        self,
        *,
        dx_um: float = 0.0,
        dy_um: float = 0.0,
        dz_um: float = 0.0,
    ) -> SimulationDomain:
        """Return a new domain with expanded dimensions.

        Args:
            dx_um: Additional length along x in um.
            dy_um: Additional half-width along y in um.
            dz_um: Additional depth along z in um.

        Returns:
            A new SimulationDomain with the specified expansions added.
        """
        return SimulationDomain(
            x_length_um=self.x_length_um + dx_um,
            y_length_um=self.y_length_um + dy_um,
            z_depth_um=self.z_depth_um + dz_um,
            spatial_resolution_um=self.spatial_resolution_um,
        )


@dataclass(frozen=True)
class MeltPoolResult:
    """Melt pool geometry and temperature extremes.

    Attributes:
        length: Melt pool length along x in metres.
        width: Melt pool full width along y in metres.
        depth: Melt pool depth along z in metres.
        peak_temperature: Maximum temperature in the domain in Kelvin.
        min_temperature: Minimum temperature in the domain in Kelvin.
    """

    length: float
    width: float
    depth: float
    peak_temperature: float
    min_temperature: float

    @property
    def length_um(self) -> float:
        """Melt pool length in micrometres."""
        return self.length / _UM_TO_M

    @property
    def width_um(self) -> float:
        """Melt pool full width in micrometres."""
        return self.width / _UM_TO_M

    @property
    def depth_um(self) -> float:
        """Melt pool depth in micrometres."""
        return self.depth / _UM_TO_M
