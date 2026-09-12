"""Exact piecewise-linear forcing convolution of finite-cylinder heat modes."""
import numpy as np
from scipy.optimize import brentq
from scipy.special import j0,j1,jn_zeros


def finite_cylinder_heat(times,radii,z,air,initial=28.,alpha=.36/(820*2600),R=.02,H=.125,h=25.,k=.36,modes=128):
    bi=h*R/k;biz=h*H/k
    lo=np.r_[0.,jn_zeros(1,modes-1)];hi=jn_zeros(0,modes)
    lam=np.array([brentq(lambda x:x*j1(x)-bi*j0(x),a,b,xtol=1e-13) for a,b in zip(lo,hi)])
    mu=np.array([brentq(lambda x:x*np.sin(x)-biz*np.cos(x),n*np.pi,n*np.pi+np.pi/2,xtol=1e-13) for n in range(modes)])
    ar=2*j1(lam)/(lam*(j0(lam)**2+j1(lam)**2))
    az=4*np.sin(mu)/(2*mu+np.sin(2*mu))
    rate=alpha*((lam[:,None]/R)**2+(mu[None,:]/H)**2)
    br=j0(lam[:,None]*np.asarray(radii)[None,:]/R)
    bz=np.cos(mu[:,None]*np.asarray(z)[None,:]/H)
    times=np.asarray(times);cuts=np.unique(np.r_[times,air[air[:,0]<=times[-1],0]])
    state=np.full_like(rate,air[0,1]-initial);out=[]
    for a,b in zip(cuts[:-1],cuts[1:]):
        slope=(np.interp(b,air[:,0],air[:,1])-np.interp(a,air[:,0],air[:,1]))/(b-a)
        decay=np.exp(-rate*(b-a));state=state*decay+slope*(-np.expm1(-rate*(b-a)))/rate
        if b in times:
            value=np.interp(b,air[:,0],air[:,1])-np.einsum("mn,mr,nz->rz",state*ar[:,None]*az[None,:],br,bz,optimize=False)
            out.append(value)
    result=np.array(out)
    if times[0]==0:result=np.concatenate((np.full((1,len(radii),len(z)),initial),result))
    assert len(result)==len(times)
    return result
