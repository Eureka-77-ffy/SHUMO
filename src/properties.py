"""Appendix parameters; thermal effective density is not dry-solid density."""
import numpy as np
from scipy.special import expi


class PropertyLaw:
    def __init__(self, config, group):
        self.group = group
        self.p = config["property_sets"][group]

    def evaluate(self, C, T_C):
        C, T_C = np.asarray(C), np.asarray(T_C)
        if np.any(C <= 0) or np.any(T_C + 273.15 <= 0):
            raise ValueError("Nonphysical accepted/trial property state: C<=0 or T_K<=0")
        p = self.p
        rho = p["rho_kg_m3"]["offset"] + p["rho_kg_m3"]["slope_C"] * C
        cp = p["cp_J_kg_K"]["offset"] + p["cp_J_kg_K"]["fraction_C_over_1_plus_C"] * C/(1+C)
        k = p["k_W_m_K"]["offset"] + p["k_W_m_K"]["fraction_C_over_1_plus_C"] * C/(1+C)
        d = p["D_m2_s"]
        D = d["prefactor"] * np.exp(-d["moisture_exponent"]/C-d["thermal_exponent_K"]/(T_C+273.15))
        return rho*cp, k, D


class Diffusivity:
    """Scalar field coefficient; Q1 D(C), or a constant validation coefficient."""
    def __init__(self, prefactor, a=0.0):
        self.prefactor, self.a = float(prefactor), float(a)

    def value(self, u):
        u = np.asarray(u)
        if self.a and np.any(u <= 0):
            raise ValueError("Nonpositive moisture state, no silent clipping")
        return np.ones_like(u, dtype=float)*self.prefactor if not self.a else self.prefactor*np.exp(-self.a/u)

    def derivative(self, u):
        return np.zeros_like(u, dtype=float) if not self.a else self.value(u)*self.a/np.asarray(u)**2

    def primitive_scaled(self, u):
        u = np.asarray(u)
        if not self.a:
            return u
        return u*np.exp(-self.a/u)+self.a*expi(-self.a/u)

    def integral(self, lo, hi):
        return self.prefactor*(self.primitive_scaled(hi)-self.primitive_scaled(lo))
