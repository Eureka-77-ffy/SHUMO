"""Moving-domain transport on a fixed reference-radius mesh.

Let s=R(t)/R0. The reference operator is multiplied by s^-2, and
reference Robin coefficients are s*h and s*hm. No Rdot advection is added:
material velocity equals grid velocity. Dry-basis state is not diluted.
"""
from pathlib import Path
import numpy as np
from openpyxl import load_workbook
from scipy.interpolate import PchipInterpolator

from src.coupled.model import CoupledFVM


class Radius:
    def __init__(self, config, kind="linear", extension="error"):
        if kind not in ("linear","pchip","fixed") or extension not in ("error","hold_last"):
            raise ValueError("Unknown radius scenario")
        self.kind,self.extension=kind,extension
        self.R0=config["geometry"]["reference_radius_m"]
        wb=load_workbook(Path(config["source_directory"])/"附件2.xlsx",read_only=True,data_only=True)
        self.data=np.asarray(list(wb.active.values)[1:],float); wb.close()
        self.times=self.data[:,0]
        self.radii=self.data[:,1]*.01
        self.curve=PchipInterpolator(self.times,self.radii,extrapolate=False)
        assert self.times[0]==0 and self.radii[0]==self.R0

    def value(self,t):
        t=np.asarray(t,dtype=float)
        if np.any(t<0):
            raise ValueError("Negative radius time")
        if self.kind=="fixed":
            return np.full_like(t,self.R0)
        if np.any(t>self.times[-1]) and self.extension=="error":
            raise ValueError("Radius observations end at 72 h; no unlabelled extrapolation")
        tq=np.minimum(t,self.times[-1])
        return np.interp(tq,self.times,self.radii) if self.kind=="linear" else self.curve(tq)


class MovingCoupled(CoupledFVM):
    def __init__(self,mesh,config,environment,radius,group="q4",perturbation=None):
        super().__init__(mesh,config,environment,group,perturbation)
        self.physical_h,self.physical_hm=self.h,self.hm
        self.radius=radius

    def geometry(self,t):
        ratio=self.radius.value(t)/self.mesh.radius
        self.h=self.physical_h*ratio
        self.hm=self.physical_hm*ratio
        return ratio

    def surface(self,t,Ti,Ci,derivatives=False):
        self.geometry(t)
        if self.A!=0:
            return super().surface(t,Ti,Ci,derivatives)
        if self.physical_hm!=0:
            raise ValueError("Zero diffusion test requires insulated water boundary")
        Ti,Ci,Ta,_=np.broadcast_arrays(Ti,Ci,*self.ambient(t))
        _,k,_,kp=self.thermal(Ci)
        hd=self.h*self.mesh.half_width
        Ts=(k*Ti+hd*Ta)/(k+hd)
        if not derivatives:
            return Ts,Ci.copy()
        return Ts,Ci.copy(),k/(k+hd),kp*(Ti-Ts)/(k+hd),np.zeros_like(Ci),np.ones_like(Ci)

    def fluxes(self,t,state,derivatives=False):
        self.geometry(t)
        if self.A!=0:
            return super().fluxes(t,state,derivatives)
        T,C=state[::2],state[1::2]
        if np.any(C<=0) or np.any(T+273.15<=0):
            raise ValueError("Nonphysical test state")
        b,k,bp,kp=self.thermal(C)
        dl,dr=(self.shaped(x,T.ndim) for x in (self.mesh.left_distances,self.mesh.right_distances))
        resistance=dl/k[:-1]+dr/k[1:]
        qt=(T[:-1]-T[1:])/resistance
        surf=self.surface(t,T[-1],C[-1],derivatives)
        qts=self.h*(surf[0]-self.ambient(t)[0])
        qc=np.zeros_like(qt); qcs=np.zeros_like(qts)
        if not derivatives:
            return b,bp,qt,qc,qts,qcs
        dqt=(1/resistance,qt/resistance*dl*kp[:-1]/k[:-1]**2,
             -1/resistance,qt/resistance*dr*kp[1:]/k[1:]**2)
        return b,bp,qt,qc,qts,qcs,dqt,(qc,qc,qc,qc),surf

    def rhs(self,t,state):
        ratio=self.geometry(t)
        return super().rhs(t,state)/ratio**2

    def jacobian(self,t,state):
        ratio=self.geometry(t)
        return super().jacobian(t,state)/ratio**2
