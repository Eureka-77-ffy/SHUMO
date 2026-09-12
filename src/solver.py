"""Piecewise adaptive stiff integration and independent flux quadrature."""
import time
import numpy as np
from scipy.integrate import solve_ivp, BDF
from numpy.polynomial.legendre import leggauss


class InitializedBDF(BDF):
    """Initialize unused work rows; leave SciPy's BDF formulas unchanged.

    SciPy 1.13.1 uses np.empty for D and initially assigns only D[0:2].
    Its first difference update reads the not-yet-used D[2], potentially
    emitting an invalid-subtract warning from arbitrary allocator contents.
    These higher differences are filled before they affect the solution.
    Explicit initialization makes warning-as-error runs deterministic.
    No warnings or nonphysical accepted states are suppressed.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.D[2:] = 0.0


def integrate_scalar(model, initial, end, sample_times, knots, method="BDF", rtol=1e-9, atol=1e-11, max_step=60., balance_order=8):
    started=time.perf_counter()
    n=model.mesh.n
    initial=np.broadcast_to(np.asarray(initial,float),(n,)).copy()
    state=np.r_[initial,0.0]
    sample_times=np.asarray(sample_times,float)
    if np.any(np.diff(sample_times)<=0) or sample_times[0]<0 or sample_times[-1]>end:
        raise ValueError("Invalid evaluation timeline")
    cuts=np.unique(np.r_[0.,np.asarray(knots)[(np.asarray(knots)>0)&(np.asarray(knots)<end)],end])
    output=np.empty((len(sample_times),n+1))
    gx,gw=leggauss(balance_order)
    accepted=[]
    quadrature_out=0.0
    max_quadrature_balance=0.0
    nfev=njev=nlu=0
    for a,b in zip(cuts[:-1],cuts[1:]):
        integrator=InitializedBDF if method=="BDF" else method
        result=solve_ivp(model.rhs,(a,b),state,method=integrator,rtol=rtol,atol=np.r_[np.full(n,atol),min(atol,1e-11)],jac=model.jacobian,dense_output=True,max_step=max_step)
        if not result.success:
            raise RuntimeError(result.message)
        if not np.isfinite(result.y).all():
            raise RuntimeError("Nonfinite accepted solution")
        state=result.y[:,-1]
        indices=np.flatnonzero((sample_times>=a)&(sample_times<=b))
        if len(indices):
            output[indices]=result.sol(sample_times[indices]).T
        steps=np.diff(result.t)
        accepted.extend(steps.tolist())
        nfev+=result.nfev; njev+=result.njev; nlu+=result.nlu
        mid=.5*(result.t[1:]+result.t[:-1]); half=.5*steps
        tq=(mid[:,None]+half[:,None]*gx).reshape(-1)
        inner=result.sol(tq)[n-1]
        rate=2*model.surface_flux(tq,inner)/(model.mesh.radius*model.capacity)
        increments=half*np.sum(rate.reshape(-1,balance_order)*gw,axis=1)
        cumulative=quadrature_out+np.cumsum(increments)
        weighted_means=np.sum(model.mesh.weights[:,None]*result.y[:n,1:],axis=0)
        residual=weighted_means-np.sum(model.mesh.weights*initial)+cumulative
        max_quadrature_balance=max(max_quadrature_balance,float(np.max(np.abs(residual))))
        quadrature_out=float(cumulative[-1])
    mean=np.sum(output[:,:n]*model.mesh.weights,axis=1)
    algebraic_balance=mean-np.sum(model.mesh.weights*initial)+output[:,-1]
    stats={"method":method,"rtol":rtol,"atol":atol,"max_step_cap_s":max_step,
           "accepted_steps":len(accepted),"accepted_max_step_s":max(accepted),
           "accepted_min_step_s":min(accepted),"nfev":nfev,"njev":njev,"nlu":nlu,
           "balance_quadrature_order":balance_order,"independent_balance_max_abs":max_quadrature_balance,
           "augmented_balance_max_abs":float(np.max(np.abs(algebraic_balance))),
           "elapsed_s":time.perf_counter()-started}
    return output,stats
