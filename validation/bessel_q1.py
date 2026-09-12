"""Independent Fourier-Bessel + exact piecewise-linear forcing convolution.

No FVM matrices, boundary recovery or numerical time integrator is used.
Eigenvalues satisfy lambda*J1(lambda)=Bi*J0(lambda).
"""
import numpy as np
from scipy.special import j0,j1,jn_zeros
from scipy.optimize import brentq


def cylinder_response(times,radii,knots,ambient_values,initial,alpha,radius,bi,modes=256):
    times,radii=np.asarray(times,float),np.asarray(radii,float)
    knots,ambient_values=np.asarray(knots,float),np.asarray(ambient_values,float)
    if np.any(times<knots[0]) or np.any(times>knots[-1]) or bi<=0:
        raise ValueError("Invalid Bessel reference interval/Bi")
    lower=np.r_[0.,jn_zeros(1,modes-1)]
    upper=jn_zeros(0,modes)
    lam=np.array([brentq(lambda z:z*j1(z)-bi*j0(z),a,b,xtol=1e-13) for a,b in zip(lower,upper)])
    rate=alpha*lam**2/radius**2
    coefficient=2*j1(lam)/(lam*(j0(lam)**2+j1(lam)**2))
    modal_at_knots=np.empty((len(knots),modes))
    modal_at_knots[0]=ambient_values[0]-initial
    slopes=np.diff(ambient_values)/np.diff(knots)
    for k,dt in enumerate(np.diff(knots)):
        decay=np.exp(-rate*dt)
        modal_at_knots[k+1]=modal_at_knots[k]*decay+slopes[k]*(-np.expm1(-rate*dt))/rate
    interval=np.clip(np.searchsorted(knots,times,side="right")-1,0,len(knots)-2)
    elapsed=times-knots[interval]
    decay=np.exp(-elapsed[:,None]*rate)
    states=modal_at_knots[interval]*decay+slopes[interval,None]*(-np.expm1(-elapsed[:,None]*rate))/rate
    basis=j0(lam[:,None]*radii[None,:]/radius)
    # Explicit contraction avoids Accelerate/NumPy 2.0 spurious floating-point
    # status warnings observed on this host; no warning suppression is used.
    result=np.interp(times,knots,ambient_values)[:,None]-np.einsum("tm,mr->tr",states*coefficient,basis,optimize=False)
    # Uniform initial state is specified strongly in the interior, before any
    # incompatible Robin initial-boundary layer in constant-diffusion tests.
    result[times==0]=initial
    return result
