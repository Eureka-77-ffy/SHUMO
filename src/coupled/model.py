"""Fixed-domain, coupled annular FVM with an analytic block Jacobian.

Surface moisture resistance integrates D(C,T_face) in C, holding the local
half-cell mean temperature fixed, not claiming a global Kirchhoff transform.
Heat resistance uses k at the half-cell mean moisture. Both surface states
are solved together. Every coefficient is updated from the current state.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix

from src.properties import PropertyLaw, Diffusivity


class CoupledFVM:
    def __init__(self, mesh, config, environment, group="q23", perturbation=None):
        self.mesh, self.config, self.environment = mesh, config, environment
        self.group = group
        self.law = PropertyLaw(config, group)
        self.p = self.law.p
        self.options = dict(perturbation or {})
        allowed = {"h_multiplier", "hm_multiplier", "D_multiplier", "beta_multiplier",
                   "freeze_D_temperature", "freeze_thermal_moisture"}
        if set(self.options) - allowed:
            raise ValueError("Unknown perturbation")
        for key, value in self.options.items():
            if key.endswith("multiplier") and value <= 0:
                raise ValueError("Nonpositive multiplier")
        self.h = config["boundary"]["heat_transfer_W_m2_K"] * self.options.get("h_multiplier", 1.)
        self.hm = config["boundary"]["normalized_mass_transfer_m_s"] * self.options.get("hm_multiplier", 1.)
        self.beta = self.options.get("beta_multiplier", 1.)
        self.C0 = config["initial"]["moisture_dry_basis"]
        self.T0 = config["initial"]["temperature_C"]
        dp = self.p["D_m2_s"]
        self.A = dp["prefactor"] * self.options.get("D_multiplier", 1.)
        self.a, self.B = dp["moisture_exponent"], dp["thermal_exponent_K"]
        self.primitive = Diffusivity(1., self.a)
        self._indices()

    def _indices(self):
        n = self.mesh.n
        i = np.arange(n-1)
        columns = (2*i, 2*i+1, 2*i+2, 2*i+3)
        self.rows = np.concatenate([np.full_like(i, 0)+row for row in
                                    (2*i, 2*i+2, 2*i+1, 2*i+3) for _ in columns]
                                   + [np.array([2*n-2, 2*n-2, 2*n-1, 2*n-1]), 2*np.arange(n)])
        self.cols = np.concatenate([col for _ in range(4) for col in columns]
                                   + [np.array([2*n-2, 2*n-1, 2*n-2, 2*n-1]), 2*np.arange(n)+1])

    def ambient(self, t):
        return self.environment.value(t, "temperature"), self.beta*self.environment.value(t, "moisture")

    def thermal(self, C):
        p = self.p
        frozen = self.options.get("freeze_thermal_moisture", False)
        c = np.full_like(C, self.C0) if frozen else C
        rho = p["rho_kg_m3"]["offset"] + p["rho_kg_m3"]["slope_C"]*c
        cp = p["cp_J_kg_K"]["offset"] + p["cp_J_kg_K"]["fraction_C_over_1_plus_C"]*c/(1+c)
        k = p["k_W_m_K"]["offset"] + p["k_W_m_K"]["fraction_C_over_1_plus_C"]*c/(1+c)
        bp = (p["rho_kg_m3"]["slope_C"]*cp
              + rho*p["cp_J_kg_K"]["fraction_C_over_1_plus_C"]/(1+c)**2)
        kp = p["k_W_m_K"]["fraction_C_over_1_plus_C"]/(1+c)**2
        if frozen:
            bp, kp = np.zeros_like(c), np.zeros_like(c)
        return rho*cp, k, bp, kp

    def diffusion(self, C, T):
        frozen = self.options.get("freeze_D_temperature", False)
        temp = np.full_like(T, self.T0) if frozen else T
        D = self.A*np.exp(-self.a/C-self.B/(temp+273.15))
        DC = D*self.a/C**2
        DT = np.zeros_like(D) if frozen else D*self.B/(temp+273.15)**2
        return D, DC, DT

    def surface(self, t, Ti, Ci, derivatives=False):
        """A bracketed scalar solve after eliminating Ts(Cs); no state clipping."""
        Ti, Ci, Ta, Ce = np.broadcast_arrays(Ti, Ci, *self.ambient(t))
        if np.any(Ci <= 0) or np.any(Ce <= 0) or np.any(Ti+273.15 <= 0):
            raise ValueError("Nonphysical surface input")
        hd, md = self.h*self.mesh.half_width, self.hm*self.mesh.half_width
        frozen = self.options.get("freeze_D_temperature", False)

        def evaluate(cs):
            _, k, _, kp = self.thermal((Ci+cs)/2)
            ts = (k*Ti+hd*Ta)/(k+hd)
            tm = np.full_like(ts, self.T0) if frozen else (Ti+ts)/2
            z = md/(self.A*np.exp(-self.B/(tm+273.15)))
            L = np.zeros_like(tm) if frozen else self.B/(tm+273.15)**2
            ts_c = .5*kp*(Ti-ts)/(k+hd)
            ts_t = k/(k+hd)
            g = self.primitive.primitive_scaled(Ci)-self.primitive.primitive_scaled(cs)-z*(cs-Ce)
            gs = -np.exp(-self.a/cs)-z + .5*z*L*ts_c*(cs-Ce)
            return ts, g, gs, z, L, ts_c, ts_t

        D = self.diffusion(Ci, Ti)[0]
        cs = (D*Ci+md*Ce)/(D+md)
        lo, hi = np.minimum(Ci, Ce).copy(), np.maximum(Ci, Ce).copy()
        for _ in range(16):
            ts, g, gs, *_ = evaluate(cs)
            lo = np.where(g > 0, cs, lo)
            hi = np.where(g < 0, cs, hi)
            proposal = cs-g/gs
            proposal = np.where((proposal < lo) | (proposal > hi), (lo+hi)/2, proposal)
            delta = np.max(np.abs(proposal-cs))
            cs = proposal
            if delta < 2e-13:
                break
        ts, g, gs, z, L, ts_c, ts_t = evaluate(cs)
        bad = np.abs(g) > 3e-12
        if np.any(bad):
            # Rare fallback uses an independently bracketed scalar root.
            flat = np.asarray(cs).reshape(-1).copy()
            for j in np.flatnonzero(np.asarray(bad).reshape(-1)):
                ti, ci, ta, ce = (float(v.flat[j]) for v in (Ti, Ci, Ta, Ce))
                def residual(c):
                    k = float(self.thermal(np.asarray((ci+c)/2))[1])
                    ts_ = (k*ti+hd*ta)/(k+hd)
                    tm_ = self.T0 if frozen else (ti+ts_)/2
                    z_ = md/(self.A*np.exp(-self.B/(tm_+273.15)))
                    return float(self.primitive.integral(c, ci)-z_*(c-ce))
                flat[j] = brentq(residual, min(ci, ce), max(ci, ce), xtol=5e-14)
            cs = flat.reshape(Ci.shape)
            ts, g, gs, z, L, ts_c, ts_t = evaluate(cs)
        if np.max(np.abs(g)) > 1e-10:
            raise RuntimeError("Surface nonlinear solve failed")
        if not derivatives:
            return ts, cs
        gi_t = .5*z*L*(1+ts_t)*(cs-Ce)
        gi_c = np.exp(-self.a/Ci)+.5*z*L*ts_c*(cs-Ce)
        cs_t, cs_c = -gi_t/gs, -gi_c/gs
        return ts, cs, ts_t+ts_c*cs_t, ts_c*(1+cs_c), cs_t, cs_c

    @staticmethod
    def shaped(vector, ndim):
        return vector.reshape((-1,)+(1,)*(ndim-1))

    def fluxes(self, t, state, derivatives=False):
        T, C = state[0::2], state[1::2]
        if np.any(C <= 0) or np.any(T+273.15 <= 0):
            raise ValueError("Nonphysical trial/accepted state; no clipping")
        b, k, bp, kp = self.thermal(C)
        D, DC, DT = self.diffusion(C, T)
        m = self.mesh
        dl, dr = (self.shaped(x, T.ndim) for x in (m.left_distances, m.right_distances))
        rk, rd = dl/k[:-1]+dr/k[1:], dl/D[:-1]+dr/D[1:]
        qt, qc = (T[:-1]-T[1:])/rk, (C[:-1]-C[1:])/rd
        surf = self.surface(t, T[-1], C[-1], derivatives)
        ta, ce = self.ambient(t)
        qts, qcs = self.h*(surf[0]-ta), self.hm*(surf[1]-ce)
        if not derivatives:
            return b, bp, qt, qc, qts, qcs
        dqt = (1/rk, qt/rk*dl*kp[:-1]/k[:-1]**2,
               -1/rk, qt/rk*dr*kp[1:]/k[1:]**2)
        dqc = (qc/rd*dl*DT[:-1]/D[:-1]**2,
               1/rd+qc/rd*dl*DC[:-1]/D[:-1]**2,
               qc/rd*dr*DT[1:]/D[1:]**2,
               -1/rd+qc/rd*dr*DC[1:]/D[1:]**2)
        return b, bp, qt, qc, qts, qcs, dqt, dqc, surf

    def rhs(self, t, state):
        b, _, qt, qc, qts, qcs = self.fluxes(t, state)
        area = self.shaped(self.mesh.faces[1:-1], qt.ndim)
        volume = self.shaped(self.mesh.volumes, qt.ndim)
        out = np.empty_like(state)
        for index, flux, surface, capacity in ((0, qt, qts, b), (1, qc, qcs, 1.)):
            total = np.concatenate((np.zeros_like(flux[:1]), area*flux,
                                    np.asarray(self.mesh.radius*surface).reshape((1,)+flux.shape[1:])), axis=0)
            out[index::2] = -np.diff(total, axis=0)/(volume*capacity)
        return out

    def jacobian(self, t, state):
        b, bp, qt, qc, qts, qcs, dqt, dqc, surf = self.fluxes(t, state, True)
        m = self.mesh
        a = m.faces[1:-1]
        values = ([ -a*d/(m.volumes[:-1]*b[:-1]) for d in dqt]
                  + [a*d/(m.volumes[1:]*b[1:]) for d in dqt]
                  + [-a*d/m.volumes[:-1] for d in dqc]
                  + [a*d/m.volumes[1:] for d in dqc])
        values.append(-m.radius/m.volumes[-1]*np.array([
            self.h*surf[2]/b[-1], self.h*surf[3]/b[-1], self.hm*surf[4], self.hm*surf[5]]))
        total = np.r_[0., a*qt, m.radius*qts]
        tdot = -np.diff(total)/(m.volumes*b)
        values.append(-tdot*bp/b)
        return coo_matrix((np.concatenate(values), (self.rows, self.cols)),
                          shape=(2*m.n, 2*m.n)).tocsc()
