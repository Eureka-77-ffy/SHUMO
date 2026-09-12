#!/usr/bin/env python3
"""Tests of conservation, analytic Jacobian, geometry and Robin reconstruction."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy.optimize._numdiff import approx_derivative
from src.fvm import RadialMesh,ScalarFVM
from src.properties import Diffusivity,PropertyLaw
from src.boundary import Environment,surface_value
from src.data_io import ROOT,load_config,save_json
from src.solver import integrate_scalar
from src.reconstruction import reconstruct
from validation.bessel_q1 import cylinder_response


def main():
    cfg=load_config();env=Environment(cfg);checks={}
    for group in ["q1","q23","q4"]:
        b,k,D=PropertyLaw(cfg,group).evaluate(2.55,28)
        assert min(b,k,D)>0
    checks["all_appendix_property_sets_positive_at_initial"]=True
    for grading in [1.,2.]:
        mesh=RadialMesh(24,.02,grading)
        assert abs(mesh.weights.sum()-1)<1e-14
        for scheme in ["harmonic","kirchhoff"]:
            law=Diffusivity(7e-9,.89)
            model=ScalarFVM(mesh,law,8e-7,lambda t:.03,flux_scheme=scheme)
            y=np.r_[2.55-.8*(mesh.centres/.02)**2,0.]
            analytic=model.jacobian(100.,y).toarray()
            finite=approx_derivative(lambda z:model.rhs(100.,z),y,method="3-point",rel_step=2e-5)
            error=np.max(np.abs(analytic-finite))/np.max(np.abs(analytic))
            assert error<1e-6,(grading,scheme,error)
            dy=model.rhs(100.,y)
            assert abs(mesh.weights@dy[:-1]+dy[-1])<1e-12
            checks[f"jacobian_g{grading}_{scheme}_relative_error"]=float(error)
        equilibrium=ScalarFVM(mesh,Diffusivity(7e-9,.89),8e-7,lambda t:2.55)
        assert np.max(np.abs(equilibrium.rhs(0,np.r_[np.full(24,2.55),0.])))==0
        insulated=ScalarFVM(mesh,Diffusivity(.36),0,lambda t:50.,820*2600)
        state=np.r_[28+2*(mesh.centres/.02)**2,0.]
        assert abs(mesh.weights@insulated.rhs(0,state)[:-1])<1e-12
    checks["constant_state_and_zero_boundary_flow"]=True
    # Quadratic radial profile has Laplacian 4 in every interior uniform cell.
    mesh=RadialMesh(40,.02)
    model=ScalarFVM(mesh,Diffusivity(1),0,lambda t:0.)
    derivative=model.rhs(0,np.r_[mesh.centres**2,0.])
    assert np.max(np.abs(derivative[:-2]-4))<1e-10
    recovered=reconstruct(model,[0],(1+mesh.centres**2)[None,:],[0])
    assert abs(recovered[0,0]-1)<1e-14
    checks["quadratic_axis_limit_and_axis_reconstruction"]=True
    for radius in [.02,.016,.01198]:
        mesh=RadialMesh(20,radius)
        off=ScalarFVM(mesh,Diffusivity(0),0,lambda t:0.)
        state=np.r_[1+(mesh.centres/radius)**2,0.]
        assert np.max(np.abs(off.rhs(0,state)))==0
    checks["no_transport_reference_profile_unchanged_for_three_radii"]=True
    ci=np.array([2.55,2.,.8]);ce=np.array([.02,.05,.03]);delta=2e-5;law=Diffusivity(7e-9,.89)
    cs=surface_value(ci,ce,delta,law,8e-7)
    residual=law.integral(cs,ci)/delta-8e-7*(cs-ce)
    assert np.max(np.abs(residual))<1e-14
    assert np.all((cs>ce)&(cs<ci))
    checks["true_surface_robin_max_flux_residual"]=float(np.max(np.abs(residual)))
    # Constant-D water transfer provides an independent mass-field benchmark.
    mesh=RadialMesh(320,.02,2.)
    model=ScalarFVM(mesh,Diffusivity(5e-9),8e-7,lambda t:.03)
    ts=np.array([0.,1.,10.,100.,600.])
    values,stats=integrate_scalar(model,2.55,600,ts,[0,600],rtol=1e-10,atol=1e-12)
    sampled=reconstruct(model,ts,values[:,:320],np.linspace(0,.02,21),uniform_initial=2.55)
    exact=cylinder_response(ts,np.linspace(0,.02,21),[0.,600.],[.03,.03],2.55,5e-9,.02,8e-7*.02/5e-9,512)
    error=float(np.max(np.abs(sampled-exact)))
    assert error<2e-4,error
    checks["constant_D_water_Bessel_max_error"]=error
    checks["constant_D_water_independent_balance"]=stats["independent_balance_max_abs"]
    save_json(ROOT/"reports/kernel_checks.json",checks)
    print(checks,flush=True)


if __name__=="__main__":
    main()
