"""Conservative annular scalar transport with analytic tridiagonal Jacobian.

    capacity * u_t = (r*a(u)*u_r)_r/r.
Only the scalar/decoupled Q1 solver is validated in this milestone.
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import diags

from .boundary import surface_value


@dataclass
class RadialMesh:
    n: int
    radius: float
    grading: float = 1.0

    def __post_init__(self):
        if self.n < 3 or self.radius <= 0 or self.grading < 1:
            raise ValueError("Invalid radial mesh")
        uniform = np.linspace(0, 1, self.n+1)
        self.faces = self.radius*(1-(1-uniform)**self.grading)
        self.centres = .5*(self.faces[1:]+self.faces[:-1])
        self.volumes = .5*np.diff(self.faces**2)  # per unit 2*pi*L
        self.weights = 2*self.volumes/self.radius**2
        self.left_distances = self.faces[1:-1]-self.centres[:-1]
        self.right_distances = self.centres[1:]-self.faces[1:-1]
        self.separations = np.diff(self.centres)
        self.half_width = self.faces[-1]-self.centres[-1]


class ScalarFVM:
    def __init__(self, mesh, coefficient, transfer, ambient, capacity=1.0, flux_scheme="harmonic"):
        self.mesh, self.coefficient = mesh, coefficient
        self.transfer, self.ambient, self.capacity = transfer, ambient, capacity
        self.flux_scheme = flux_scheme
        if flux_scheme not in ("harmonic", "kirchhoff"):
            raise ValueError("Unknown face flux")
        if transfer < 0 or capacity <= 0 or coefficient.prefactor < 0:
            raise ValueError("Invalid transport coefficient/capacity")

    def boundary(self, t, inner):
        return surface_value(inner, self.ambient(t), self.mesh.half_width, self.coefficient, self.transfer)

    def surface_flux(self, t, inner):
        return self.transfer*(self.boundary(t, inner)-self.ambient(t))

    def face_flux_and_derivatives(self, t, u):
        m, law = self.mesh, self.coefficient
        if law.prefactor == 0:
            if self.transfer != 0:
                raise ValueError("Zero diffusion test requires insulated surface")
            zeros = np.zeros(m.n-1)
            return zeros, zeros.copy(), zeros.copy(), 0.0, 0.0
        D = law.value(u)
        if self.flux_scheme == "kirchhoff":
            flux = -law.integral(u[:-1], u[1:])/m.separations
            left, right = D[:-1]/m.separations, -D[1:]/m.separations
        else:
            resistance = m.left_distances/D[:-1]+m.right_distances/D[1:]
            flux = (u[:-1]-u[1:])/resistance
            Dprime = law.derivative(u)
            left = 1/resistance + flux/resistance*m.left_distances*Dprime[:-1]/D[:-1]**2
            right = -1/resistance + flux/resistance*m.right_distances*Dprime[1:]/D[1:]**2
        surface = self.boundary(t,u[-1])
        q_s = self.transfer*(surface-self.ambient(t))
        derivative_s = self.transfer*D[-1]/(law.value(surface)+self.transfer*m.half_width)
        return flux, left, right, q_s, derivative_s

    def rhs(self, t, state):
        n, m = self.mesh.n, self.mesh
        flux, _, _, qs, _ = self.face_flux_and_derivatives(t, state[:n])
        face_total = np.concatenate(([0.0], m.faces[1:-1]*flux, [m.radius*qs]))
        du = -np.diff(face_total)/(m.volumes*self.capacity)
        cumulative_outward = 2*qs/(m.radius*self.capacity)
        return np.r_[du, cumulative_outward]

    def jacobian(self, t, state):
        m, n = self.mesh, self.mesh.n
        _, left, right, _, ds = self.face_flux_and_derivatives(t, state[:n])
        area = m.faces[1:-1]
        scale = m.volumes*self.capacity
        diagonal = np.zeros(n+1)
        diagonal[:-2] -= area*left/scale[:-1]
        diagonal[1:-1] += area*right/scale[1:]
        diagonal[n-1] -= m.radius*ds/scale[-1]
        lower = np.r_[area*left/scale[1:], 2*ds/(m.radius*self.capacity)]
        upper = np.r_[-area*right/scale[:-1], 0.0]
        return diags([lower, diagonal, upper], [-1,0,1], format="csc")
