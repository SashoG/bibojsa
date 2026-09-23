"""Hard-coded Sellmeier material data for validation and target crystals.

Each material module exports a dict of dispersion coefficients keyed by
dielectric-frame axis ('x','y','z') plus metadata (validity range,
literature source), giving a uniform interface for `bibojsa.biaxial`
regardless of whether the crystal is biaxial (BIBO, KTP, ...) or uniaxial
(KDP, ...; represented as the biaxial degenerate case nx=ny or ny=nz).
"""
from .bibo import BIBO
from .kdp import KDP
from .ktp import KTP

__all__ = ["BIBO", "KDP", "KTP"]
