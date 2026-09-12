"""Fresh-start moving Q4 with true physical-coordinate output and balances."""
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from src.events.integration import global_maximum
from src.solver import InitializedBDF


def reconstruct(model,times,state):
    times=np.asarray(times)
    T,C=state[::2].T,state[1::2].T
    R=model.radius.value(times)
    surfaces=model.surface(times,T[:,-1],C[:,-1])
    R0=model.mesh.radius
    centres=model.mesh.centres/R0
    x0,x1=centres[:2]
    points=np.r_[0.,centres,1.]
    radii=np.arange(21)*.001
    values=[]; extras={"radius_current_m":R}
    for name,cells,surface,initial in zip(("temperature","moisture"),(T,C),surfaces,(model.T0,model.C0)):
        axis=(x1*x1*cells[:,0]-x0*x0*cells[:,1])/(x1*x1-x0*x0)
        nodes=np.column_stack((axis,cells,surface))
        nodes[times==0]=initial
        array=np.full((len(times),22),np.nan)
        for i in range(len(times)):
            valid=radii<=R[i]
            xi=radii[valid]/R[i]
            curve=PchipInterpolator(points,nodes[i],extrapolate=False)
            sampled=curve(xi)
            near=xi<x0
            sampled[near]=nodes[i,0]+(cells[i,1]-cells[i,0])/(x1*x1-x0*x0)*xi[near]**2
            array[i,np.flatnonzero(valid)]=sampled
            array[i,-1]=nodes[i,-1]
        values.append(array)
        extras[name+"_mean"]=np.sum(cells*model.mesh.weights,axis=1)
        if name=="moisture":
            extras["moisture_global_max"]=np.max(nodes,axis=1)
            extras["moisture_argmax_xi"]=points[np.argmax(nodes,axis=1)]
            extras["moisture_radial_increase_max"]=np.maximum(0,np.max(np.diff(nodes,axis=1),axis=1))
    return values,extras


def integrate(model,limit=259200.,rtol=1e-10,atol_temperature=1e-11,
              atol_moisture=1e-13,max_step=30.,method="BDF",output_step=60.):
    started=time.perf_counter()
    n=model.mesh.n
    initial=np.tile([model.T0,model.C0],n)
    state=initial.copy()
    # Explicitly demand supported geometry for the entire requested horizon.
    model.radius.value(limit)
    cuts=np.unique(np.r_[0.,model.environment.times,model.radius.times,
                         np.arange(0.,limit,1800.),limit])
    cuts=cuts[(cuts>=0)&(cuts<=limit)]
    times,fields,extras=[],[[],[]],{}
    profile_t,profile_y=[0.],[initial.copy()]
    accepted_t,accepted_max,steps_all=[],[],[]
    gx,gw=leggauss(8)
    water_out=heat_out=water_error=heat_error=0.
    cap0=float(model.thermal(np.asarray(model.C0))[0])
    counts={"nfev":0,"njev":0,"nlu":0}
    event_time=None; event_state=None
    minC=model.C0
    atol=np.tile([atol_temperature,atol_moisture],n)
    solver=InitializedBDF if method=="BDF" else method

    def event(t,y):
        return float(global_maximum(model,t,y)-model.config["events"]["moisture_threshold"])
    event.direction,event.terminal=-1,True

    for a,b in zip(cuts[:-1],cuts[1:]):
        model.environment.right_at_join=a>=14400.
        result=solve_ivp(model.rhs,(a,b),state,method=solver,rtol=rtol,atol=atol,
                         jac=model.jacobian,dense_output=True,max_step=max_step,events=event)
        if not result.success or not np.isfinite(result.y).all():
            raise RuntimeError(result.message)
        state=result.y[:,-1].copy()
        stop=float(result.t[-1])
        minC=min(minC,float(np.min(result.y[1::2])))
        if minC<=0:
            raise RuntimeError("Nonpositive accepted moisture")
        dt=np.diff(result.t); steps_all.extend(dt.tolist())
        for k in counts: counts[k]+=getattr(result,k)
        checks=np.sort(np.r_[result.t,.5*(result.t[:-1]+result.t[1:])])
        maxima=global_maximum(model,checks,result.sol(checks))
        accepted_t.extend(checks.tolist()); accepted_max.extend(maxima.tolist())
        ts=np.arange(np.floor(a/output_step)*output_step+output_step,stop+1e-9,output_step)
        ts=ts[(ts>a)&(ts<=stop)]
        if a==0:
            ts=np.unique(np.r_[0.,[v for v in (1.,10.) if v<=stop],ts])
        if result.t_events[0].size:
            event_time=float(result.t_events[0][0]); event_state=result.y_events[0][0].copy()
            if len(ts) and abs(ts[-1]-event_time)<1e-9: ts[-1]=event_time
            else: ts=np.r_[ts,event_time]
        elif b==limit and (not len(ts) or ts[-1]!=b):
            ts=np.r_[ts,b]
        if len(ts):
            sampled,diag=reconstruct(model,ts,result.sol(ts))
            times.append(ts)
            for j in range(2): fields[j].append(sampled[j])
            for k,v in diag.items(): extras.setdefault(k,[]).append(v)
        for first in range(0,len(dt),32):
            last=min(first+32,len(dt)); widths=dt[first:last]
            mid=.5*(result.t[first+1:last+1]+result.t[first:last])
            tq=(mid[:,None]+.5*widths[:,None]*gx).ravel()
            yq=result.sol(tq); dy=model.rhs(tq,yq)
            st,sc=model.surface(tq,yq[-2],yq[-1]); ta,ce=model.ambient(tq)
            R=model.radius.value(tq)
            correction=np.sum(model.mesh.weights[:,None]*model.thermal(yq[1::2])[2]*(yq[::2]-model.T0)*dy[1::2],axis=0)
            heat_rate=2*model.physical_h/R*(st-ta)-correction
            water_rate=2*model.physical_hm/R*(sc-ce)
            hi=heat_out+np.cumsum(.5*widths*np.sum(heat_rate.reshape(-1,8)*gw,axis=1))
            wi=water_out+np.cumsum(.5*widths*np.sum(water_rate.reshape(-1,8)*gw,axis=1))
            ye=result.y[:,first+1:last+1]
            H=np.sum(model.mesh.weights[:,None]*model.thermal(ye[1::2])[0]*(ye[::2]-model.T0),axis=0)
            Cmean=np.sum(model.mesh.weights[:,None]*ye[1::2],axis=0)
            heat_error=max(heat_error,float(np.max(np.abs(H+hi))))
            water_error=max(water_error,float(np.max(np.abs(Cmean-model.C0+wi))))
            heat_out,water_out=float(hi[-1]),float(wi[-1])
        profile_t.append(stop); profile_y.append(state.copy())
        if event_time is not None: break
        if b%21600==0:
            print(f"Q4 {model.group}/{model.radius.kind}, {b/3600:.0f} h: R={model.radius.value(b)*100:.5f} cm, Cmax={maxima[-1]:.7f}",flush=True)
    cache={"time_s":np.concatenate(times),"radius_fixed_m":np.arange(21)*.001,
        "temperature":np.concatenate(fields[0]),"moisture":np.concatenate(fields[1]),
        "profile_time_s":np.asarray(profile_t),"profile_state":np.asarray(profile_y),
        "final_state":state,"cell_centres_reference_m":model.mesh.centres,
        "reference_faces_m":model.mesh.faces,"weights":model.mesh.weights,
        "accepted_check_time_s":np.asarray(accepted_t),"accepted_check_global_max":np.asarray(accepted_max)}
    cache.update({k:np.concatenate(v) for k,v in extras.items()})
    stats={**counts,"accepted_steps":len(steps_all),"accepted_max_step_s":max(steps_all),
        "accepted_min_step_s":min(steps_all),"minimum_accepted_moisture":minC,
        "water_balance_max_abs":water_error,"water_balance_relative":water_error/model.C0,
        "heat_corrected_balance_max_abs_J_m3":heat_error,"heat_corrected_balance_relative":heat_error/(cap0*22),
        "heat_balance_definition":"normalized mean[b*(T-T0)] + integral(2*qT/R - mean[bprime*(T-T0)*Cdot]); no duplicate geometry term",
        "event_found":event_time is not None,"event_time_s":event_time,
        "end_time_s":stop,"radius_end_m":float(model.radius.value(stop)),
        "global_max_at_end":float(cache["moisture_global_max"][-1]),
        "domain_supported_through_s":float(model.radius.times[-1]),
        "radius_extension":model.radius.extension,"global_max_increase_at_checks":float(max(0,np.max(np.diff(accepted_max))))}
    if event_time is not None:
        after=event_time+1.
        if model.radius.kind!="fixed" and model.radius.extension=="error" and after>model.radius.times[-1]:
            raise RuntimeError("Insufficient observed radius for post-event verification")
        post=solve_ivp(model.rhs,(event_time,after),event_state,method=solver,rtol=rtol,
            atol=atol,jac=model.jacobian,dense_output=True,max_step=min(1.,max_step))
        if not post.success: raise RuntimeError("Post-event solve failed")
        before=max(result.t[0],event_time-1.)
        integer=float(np.floor(event_time)+1)
        values=[event(before,result.sol(before)),event(event_time,event_state),
                event(integer,post.sol(integer)),event(after,post.y[:,-1])]
        assert values[0]>0 and values[2]<0 and values[3]<0
        stats.update({"post_event_time_s":[before,event_time,integer,after],"post_event_g":values,
            "event_slope_per_s":(values[-1]-values[0])/(after-before),"first_strict_integer_second":integer})
        cache["event_state"]=event_state; cache["post_event_state"]=post.y[:,-1]
    stats["elapsed_s"]=time.perf_counter()-started
    return cache,stats
