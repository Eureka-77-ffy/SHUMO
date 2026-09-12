"""Geometry audit kernels: degeneration, Jacobian, balances, exact 2D heat."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from src.data_io import ROOT,load_config,digest,save_json,verify_sources
from src.fvm import RadialMesh
from src.moving.model import MovingCoupled
from src.boundary import Environment
from validation.endface_model import AxisymmetricFVM
from validation.endface_run import fingerprint
from validation.endface_heat_reference import finite_cylinder_heat
from validation.bessel_q1 import cylinder_response


def main():
    verify_sources();cfg=load_config();tests=[]
    for group in ("q1","q23","q4"):
        m=AxisymmetricFVM(cfg,group,12,9,1.);t=10000.
        rr,zz=np.meshgrid(m.rc,m.zc/m.H);y=np.empty(2*m.n)
        y[::2]=(45+rr+zz).ravel();y[1::2]=(1.5-.5*rr**2-.1*zz**2).ravel()
        v=np.random.default_rng(20).normal(size=y.size);j=m.jacobian(t,y)
        fd=(m.rhs(t,y+1e-5*v)-m.rhs(t,y-1e-5*v))/2e-5
        jac=float(np.linalg.norm(fd-j@v)/np.linalg.norm(fd));assert jac<1e-7
        dy,rate=m.evaluate(t,y);balance=abs(float(m.weights@dy[1::2])+rate);assert balance<1e-12
        m.end_multiplier=0.;y[::2]=np.tile(45+m.rc,9);y[1::2]=np.tile(1.5-.5*m.rc**2,9)
        ref=MovingCoupled(RadialMesh(12,.02,2),cfg,m.environment,m.radius,group)
        delta=float(np.max(abs(m.rhs(t,y).reshape(9,24)-ref.rhs(t,y[:24]))));assert delta<1e-10
        tests.append({"group":group,"jacobian_relative_error":jac,"two_boundary_water_balance_rate_error":balance,"insulated_end_1D_degeneration_error":delta})
    # Variable coefficients and both directions: evaluate an analytic smooth
    # manufactured spatial operator away from the physical boundaries. This
    # tests the 1/R² vs 1/H² factors, cross derivatives and cylindrical 1/xi.
    manufactured=[]
    for n in (16,32,64,128):
        m=AxisymmetricFVM(cfg,"q4",n,n,1.);t=10000.;R=float(m.radius.value(t));H=m.H
        x,z=np.meshgrid(m.rc,m.zc/H);T=40+x*x+2*z*z;C=2.5-.1*x*x-.15*z*z
        y=np.column_stack((T.ravel(),C.ravel())).ravel();dy=m.rhs(t,y).reshape(n,n,2)
        b,k,_,kp=m.helper.thermal(C);D,DC,DT=m.helper.diffusion(C,T)
        dot=(-.2*x/R)*(2*x/R)+(-.3*z/H)*(4*z/H)
        normC=(.2*x/R)**2+(.3*z/H)**2
        exactT=(k*(4/R**2+4/H**2)+kp*dot)/b
        exactC=D*(-.4/R**2-.3/H**2)+DC*normC+DT*dot
        # Fixed interior physical subregion avoids moving the error norm onto
        # the imposed Robin mismatch of this pointwise manufactured test.
        mask=(x>.15)&(x<.75)&(z>.15)&(z<.75)
        manufactured.append({"n":n,"max_T_operator_error":float(np.max(abs(dy[:,:,0][mask]-exactT[mask]))),
                              "max_C_operator_error":float(np.max(abs(dy[:,:,1][mask]-exactC[mask])))})
    for key in ("max_T_operator_error","max_C_operator_error"):
        assert manufactured[-1][key]<manufactured[-2][key]/3
    env=Environment(cfg);ts=np.unique(np.r_[0.,100.,np.arange(60.,1801.,60.)]);rr=np.arange(21)*.001;zs=np.array([0.,.0625,.1125,.125])
    h256=finite_cylinder_heat(ts,rr,zs,env.data,modes=256)
    h512=finite_cylinder_heat(ts,rr,zs,env.data,modes=512)
    infinite=cylinder_response(ts,rr,env.times,env.data[:,1],28.,.36/(820*2600),.02,25*.02/.36,512)
    modal=float(np.max(abs(h512-h256)));assert modal<5e-6
    heat={"time_s":ts.tolist(),"radius_m":rr.tolist(),"z_m":zs.tolist(),"finite_temperature_C":h512.tolist(),
          "modes_per_direction":512,"256_to_512_max_change_C":modal,
          "midplane_256_to_512_change_C":float(np.max(abs(h512[:,:,0]-h256[:,:,0]))),
          "midplane_finite_minus_infinite_max_C":float(np.max(abs(h512[:,:,0]-infinite))),
          "centre_end_temperature_at_1800_C":float(h512[-1,0,-1]),"centre_mid_temperature_at_1800_C":float(h512[-1,0,0])}
    cache_comparison=[]
    for nr,nz in ((32,32),(64,32),(64,64),(128,64)):
        p=ROOT/f"results/cache/endface/q1_r{nr}_z{nz}_e1_tol1e-08_dt120.npz"
        if not p.exists():continue
        meta=json.loads(p.with_suffix(".json").read_text());assert meta["fingerprint"]==fingerprint() and meta["cache_sha256"]==digest(p)
        a=dict(np.load(p));ta=a["time_s"]
        exact=finite_cylinder_heat(ta,rr,[0.,.125],env.data,modes=512)
        mid=float(np.max(abs(a["temperature_mid"][:,:21]-exact[:,:,0])))
        delta=abs(a["temperature_end"][:,:20]-exact[:,:20,1])
        cache_comparison.append({"nr":nr,"nz":nz,"midplane_heat_max_difference_C":mid,"end_heat_max_difference_C_excluding_corner":float(np.max(delta))})
    report={"status":"passed_endface_geometry_kernels_and_exact_heat_reference","fingerprint":fingerprint(),"script_sha256":digest(__file__),
            "heat_reference_script_sha256":digest(ROOT/"validation/endface_heat_reference.py"),"kernel_cases":tests,
            "manufactured_interior_operator_checks":manufactured,"not_full_manufactured_initial_boundary_value_solve":True,
            "exact_linear_heat_reference":heat,"Q1_FVM_vs_exact":cache_comparison}
    save_json(ROOT/"reports/endface_kernel_checks.json",report)
    print(json.dumps({"manufactured":manufactured,"heat":{k:v for k,v in heat.items() if not isinstance(v,list)},"FVM":cache_comparison}),flush=True)


if __name__=="__main__":main()
