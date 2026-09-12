"""Conditional 2D axisymmetric end-face audit; never replaces official outputs.

Quarter meridian: xi=r/R(t), 0<=z<=L/2. Material v_r=Rdot*xi,
v_z=0; reference motion cancels the radial material velocity only.
Both radial and axial fluxes are retained. Shares the validated 1D local
constitutive/half-cell boundary closure, so this is a geometry-model
comparison, NOT a new independent spatial validation of that closure.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.interpolate import PchipInterpolator

from src.coupled.model import CoupledFVM
from src.fvm import RadialMesh
from src.moving.model import Radius
from src.events.integration import ScenarioEnvironment


class AxisymmetricFVM:
    def __init__(self, cfg, group, nr, nz, end_multiplier=1.):
        self.cfg,self.group,self.nr,self.nz=cfg,group,nr,nz
        self.environment=ScenarioEnvironment(cfg)
        self.radius=Radius(cfg,"linear" if group=="q4" else "fixed")
        self.H=cfg["geometry"]["length_m"]/2
        self.rf=1-(1-np.linspace(0,1,nr+1))**2
        self.zf=self.H*(1-(1-np.linspace(0,1,nz+1))**2)
        self.rc=(self.rf[1:]+self.rf[:-1])/2
        self.zc=(self.zf[1:]+self.zf[:-1])/2
        self.rvol=np.diff(self.rf**2)/2
        self.dz=np.diff(self.zf)
        self.n=nr*nz
        self.weights=(self.dz[:,None]*self.rvol[None,:]).ravel()*2/self.H
        self.helper=CoupledFVM(RadialMesh(3,.02),cfg,self.environment,group)
        self.h=cfg["boundary"]["heat_transfer_W_m2_K"]
        self.hm=cfg["boundary"]["normalized_mass_transfer_m_s"]
        self.end_multiplier=float(end_multiplier)
        ids=np.arange(self.n).reshape(nz,nr)
        self.l=np.r_[ids[:,:-1].ravel(),ids[:-1,:].ravel()]
        self.r=np.r_[ids[:,1:].ravel(),ids[1:,:].ravel()]
        self.ne_rad=nz*(nr-1)
        self.dl0=np.r_[np.tile(self.rf[1:-1]-self.rc[:-1],nz),np.repeat(self.zf[1:-1]-self.zc[:-1],nr)]
        self.dr0=np.r_[np.tile(self.rc[1:]-self.rf[1:-1],nz),np.repeat(self.zc[1:]-self.zf[1:-1],nr)]
        self.area0=np.r_[(self.dz[:,None]*self.rf[None,1:-1]).ravel(),np.tile(self.rvol,nz-1)]
        self.side=ids[:,-1];self.end=ids[-1,:]
        cols=(2*self.l,2*self.l+1,2*self.r,2*self.r+1)
        rows=(2*self.l,2*self.r,2*self.l+1,2*self.r+1)
        self.jrows=[row for row in rows for _ in cols]
        self.jcols=[col for _ in rows for col in cols]
        for ids in (self.side,self.end):
            self.jrows += [2*ids,2*ids,2*ids+1,2*ids+1]
            self.jcols += [2*ids,2*ids+1,2*ids,2*ids+1]
        self.jrows += [2*np.arange(self.n)];self.jcols += [2*np.arange(self.n)+1]
        self.jrows=np.concatenate(self.jrows);self.jcols=np.concatenate(self.jcols)

    def geometry(self,t):
        R=float(self.radius.value(t))
        dist=np.r_[np.full(self.ne_rad,R),np.ones(len(self.l)-self.ne_rad)]
        area_scale=np.r_[np.full(self.ne_rad,R),np.full(len(self.l)-self.ne_rad,R*R)]
        volume=(R*R*self.dz[:,None]*self.rvol[None,:]).ravel()
        return R,self.dl0*dist,self.dr0*dist,self.area0*area_scale,volume

    def boundary(self,t,T,C,which,derivatives=False):
        m=self.helper
        if which=="side":
            factor=1.;delta=float(self.radius.value(t))*(1-self.rc[-1])
        else:
            factor=self.end_multiplier;delta=self.H-self.zc[-1]
        if factor==0:
            if not derivatives:return T.copy(),C.copy()
            return T.copy(),C.copy(),np.ones_like(T),np.zeros_like(T),np.zeros_like(T),np.ones_like(T)
        m.h=self.h*factor;m.hm=self.hm*factor;m.mesh.half_width=delta
        return m.surface(t,T,C,derivatives)

    def evaluate(self,t,y,jac=False):
        T,C=y[::2],y[1::2]
        if np.any(C<=0) or np.any(T+273.15<=0):raise ValueError("Nonphysical 2D state")
        R,dl,dr,areas,V=self.geometry(t)
        b,k,bp,kp=self.helper.thermal(C)
        D,DC,DT=self.helper.diffusion(C,T)
        l,r=self.l,self.r
        rk=dl/k[l]+dr/k[r];rd=dl/D[l]+dr/D[r]
        qt=(T[l]-T[r])/rk;qc=(C[l]-C[r])/rd
        def divergence(q):
            return np.bincount(r,weights=areas*q,minlength=self.n)-np.bincount(l,weights=areas*q,minlength=self.n)
        ht,mc=divergence(qt),divergence(qc)
        ta,ce=self.helper.ambient(t);bound=[]
        outward_water=0.
        for ids,which,area,factor in ((self.side,"side",R*self.dz,1.),(self.end,"end",R*R*self.rvol,self.end_multiplier)):
            surf=self.boundary(t,T[ids],C[ids],which,jac)
            qT=self.h*factor*(surf[0]-ta);qC=self.hm*factor*(surf[1]-ce)
            ht[ids]-=area*qT;mc[ids]-=area*qC
            outward_water+=float(np.sum(area*qC)/(.5*R*R*self.H))
            bound.append((ids,area,factor,surf))
        out=np.empty_like(y);out[::2]=ht/(V*b);out[1::2]=mc/V
        if not jac:return out,outward_water
        dqt=(1/rk,qt/rk*dl*kp[l]/k[l]**2,-1/rk,qt/rk*dr*kp[r]/k[r]**2)
        dqc=(qc/rd*dl*DT[l]/D[l]**2,1/rd+qc/rd*dl*DC[l]/D[l]**2,
             qc/rd*dr*DT[r]/D[r]**2,-1/rd+qc/rd*dr*DC[r]/D[r]**2)
        vals=[-areas*d/(V[l]*b[l]) for d in dqt]+[areas*d/(V[r]*b[r]) for d in dqt]
        vals += [-areas*d/V[l] for d in dqc]+[areas*d/V[r] for d in dqc]
        for ids,area,factor,s in bound:
            vals += [-area*self.h*factor*s[2]/(V[ids]*b[ids]),-area*self.h*factor*s[3]/(V[ids]*b[ids]),
                     -area*self.hm*factor*s[4]/V[ids],-area*self.hm*factor*s[5]/V[ids]]
        vals += [-out[::2]*bp/b]
        return coo_matrix((np.concatenate(vals),(self.jrows,self.jcols)),shape=(2*self.n,2*self.n)).tocsc()

    def rhs(self,t,y):return self.evaluate(t,y)[0]
    def jacobian(self,t,y):return self.evaluate(t,y,True)

    @staticmethod
    def axis_value(field,x,axis):
        return (x[1]**2*np.take(field,0,axis=axis)-x[0]**2*np.take(field,1,axis=axis))/(x[1]**2-x[0]**2)

    def full_nodes(self,t,y):
        T=y[::2].reshape(self.nz,self.nr);C=y[1::2].reshape(self.nz,self.nr)
        side=self.boundary(t,T[:,-1],C[:,-1],"side")
        end=self.boundary(t,T[-1,:],C[-1,:],"end")
        fields=[]
        for f,s,e in zip((T,C),side,end):
            # PCHIP interior reconstruction is bounded by these nodes. Axes
            # use the even-quadratic limit. A corner belongs to both boundaries;
            # bracket its value by the two recovered adjacent boundary values.
            core=np.column_stack((self.axis_value(f,self.rc,1),f,s))
            middle=self.axis_value(core,self.zc,0)
            endline=np.r_[self.axis_value(e,self.rc,0),e,max(e[-1],s[-1])]
            fields.append(np.vstack((middle,core,endline)))
        return fields

    def maximum(self,t,y):
        # Also include both corner traces: no arbitrary corner interpolation
        # can hide a larger observed boundary value in the event search.
        nodes=self.full_nodes(t,y)[1]
        return float(np.max(nodes))

    def sample(self,t,y):
        nodes=self.full_nodes(t,y);R=float(self.radius.value(t))
        x=np.r_[0.,self.rc,1.];rr=np.arange(21)*.001
        fields=[]
        for f in nodes:
            rows=[]
            for profile in (f[0],f[-1]):
                value=np.full(22,np.nan);inside=rr<=R
                value[np.flatnonzero(inside)]=PchipInterpolator(x,profile)(rr[inside]/R)
                value[-1]=profile[-1];rows.append(value)
            fields.append(rows)
        # The corner needs a separate two-boundary reconstruction. Its upper
        # moisture bracket is enough for the maximum test, but do not report
        # a guessed corner temperature/moisture as a measured surface value.
        fields[0][1][-1]=np.nan;fields[1][1][-1]=np.nan
        if rr[-1]==R:fields[0][1][20]=np.nan;fields[1][1][20]=np.nan
        return {"temperature_mid":fields[0][0],"temperature_end":fields[0][1],
                "moisture_mid":fields[1][0],"moisture_end":fields[1][1],
                "mean_C":float(self.weights@y[1::2]),"max_C":self.maximum(t,y),
                "argmax_node":np.unravel_index(np.argmax(nodes[1]),nodes[1].shape),"radius_m":R}
