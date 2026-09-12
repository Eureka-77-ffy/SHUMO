"""Measured environment and true surface values from half-cell + Robin closure."""
from pathlib import Path

import numpy as np
from openpyxl import load_workbook
from scipy.optimize import brentq


class Environment:
    def __init__(self, config):
        path = Path(config["source_directory"]) / "附件1.xlsx"
        wb = load_workbook(path, read_only=True, data_only=True)
        self.data = np.array(list(wb.active.values)[1:], dtype=float)
        wb.close()
        self.times = self.data[:, 0]
        boundary = config["boundary"]
        a, b = boundary["mean_window_s"]
        self.platform = self.data[(self.times >= a) & (self.times <= b), 1:].mean(axis=0)
        self.beta = boundary["beta"]

    def value(self, t, field):
        t = np.asarray(t)
        if np.any(t < 0):
            raise ValueError("Environment time precedes initial state")
        col = 1 if field == "temperature" else 2
        value = np.interp(t, self.times, self.data[:, col])
        value = np.where(t > self.times[-1], self.platform[col-1], value)
        return value if col == 1 else value*self.beta


def surface_value(inner, ambient, half_width, coefficient, transfer):
    """Integral of D across half-cell equals Robin flux times half-width.

    This integral-resistance closure is also used with harmonic interior fluxes.
    It differs from assuming the last cell centre is the physical surface.
    """
    inner, ambient = np.broadcast_arrays(np.asarray(inner, float), np.asarray(ambient, float))
    if transfer == 0:
        return inner.copy()
    if coefficient.a == 0:
        z = transfer*half_width/coefficient.prefactor
        return (inner + z*ambient)/(1+z)
    if np.any(inner <= 0) or np.any(ambient <= 0):
        raise ValueError("Moisture surface bracket must remain positive")
    z = transfer*half_width/coefficient.prefactor
    lo, hi = np.minimum(inner, ambient).copy(), np.maximum(inner, ambient).copy()
    initial_D = coefficient.value(inner)/coefficient.prefactor
    surface = (initial_D*inner+z*ambient)/(initial_D+z)
    primitive_inner = coefficient.primitive_scaled(inner)
    for _ in range(10):
        residual = primitive_inner-coefficient.primitive_scaled(surface)-z*(surface-ambient)
        deriv = -coefficient.value(surface)/coefficient.prefactor-z
        lo = np.where(residual > 0, surface, lo)
        hi = np.where(residual < 0, surface, hi)
        proposal = surface-residual/deriv
        # Allow exact root/endpoints; replace only genuinely out-of-bracket iterates.
        proposal = np.where((proposal < lo) | (proposal > hi), .5*(lo+hi), proposal)
        if np.max(np.abs(proposal-surface)) < 2e-13:
            surface = proposal
            break
        surface = proposal
    residual = primitive_inner-coefficient.primitive_scaled(surface)-z*(surface-ambient)
    if np.any(np.abs(residual) > 3e-12):
        flat = surface.reshape(-1)
        for i in np.flatnonzero(np.abs(residual.reshape(-1)) > 3e-12):
            ci, ce = inner.flat[i], ambient.flat[i]
            flat[i] = brentq(lambda s: coefficient.primitive_scaled(ci)-coefficient.primitive_scaled(s)-z*(s-ce), min(ci,ce), max(ci,ce), xtol=5e-14)
        surface = flat.reshape(inner.shape)
    return surface
