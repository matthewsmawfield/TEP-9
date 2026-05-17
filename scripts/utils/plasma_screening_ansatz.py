"""
Phenomenological plasma screening ansätze for Step 017 / Step 041 stress tests.

All forms are explicit placeholders pending a first-principles scalar–plasma
coupling derived from the TEP action.
"""

from __future__ import annotations

import math
from typing import Literal

PlasmaAnsatzId = Literal["exponential", "powerlaw", "none", "saturating"]


def screening_exponential(plasma_density_cm3: float, n_ref_cm3: float) -> float:
    """S = exp(-n_e / n_ref), n_ref > 0."""
    if n_ref_cm3 <= 0:
        raise ValueError("n_ref_cm3 must be positive")
    if plasma_density_cm3 <= 0:
        return 1.0
    return float(math.exp(-plasma_density_cm3 / n_ref_cm3))


def screening_powerlaw(
    plasma_density_cm3: float, n_ref_cm3: float, alpha: float = 0.3
) -> float:
    """S = (1 + n_e / n_ref)^(-alpha)."""
    if n_ref_cm3 <= 0:
        raise ValueError("n_ref_cm3 must be positive")
    if alpha < 0:
        raise ValueError("alpha must be non-negative")
    if plasma_density_cm3 <= 0:
        return 1.0
    return float((1.0 + plasma_density_cm3 / n_ref_cm3) ** (-alpha))


def screening_saturating(plasma_density_cm3: float, n_ref_cm3: float) -> float:
    """S = 1 / (1 + (n_e / n_ref)^2); bounded, smooth saturating attenuation."""
    if n_ref_cm3 <= 0:
        raise ValueError("n_ref_cm3 must be positive")
    if plasma_density_cm3 <= 0:
        return 1.0
    x = plasma_density_cm3 / n_ref_cm3
    return float(1.0 / (1.0 + x * x))


def plasma_screening_for_ansatz(
    plasma_density_cm3: float,
    n_ref_cm3: float,
    ansatz: PlasmaAnsatzId,
    *,
    powerlaw_alpha: float = 0.3,
) -> float:
    if ansatz == "exponential":
        return screening_exponential(plasma_density_cm3, n_ref_cm3)
    if ansatz == "powerlaw":
        return screening_powerlaw(
            plasma_density_cm3, n_ref_cm3, alpha=powerlaw_alpha
        )
    if ansatz == "none":
        return 1.0
    if ansatz == "saturating":
        return screening_saturating(plasma_density_cm3, n_ref_cm3)
    raise ValueError(f"unknown plasma ansatz: {ansatz!r}")
