#!/usr/bin/env python3
"""Compute mesh series and independent time-boundary spectral reference."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
from scripts.question1 import run_case
from src.data_io import ROOT,load_config,save_json
from src.boundary import Environment
from validation.bessel_q1 import cylinder_response


def main():
    env=Environment(load_config());times=np.arange(1801.);radii=np.arange(21)*.001
    ref256=cylinder_response(times,radii,env.times,env.data[:,1],28,.36/(820*2600),.02,25*.02/.36,256)
    ref512=cylinder_response(times,radii,env.times,env.data[:,1],28,.36/(820*2600),.02,25*.02/.36,512)
    spectral_error=float(np.max(np.abs(ref256-ref512)))
    print('Bessel 256 vs 512 max',spectral_error,flush=True)
    result=[]
    for grading,ns in [(1.,[80,160,320,640]),(2.,[80,160,320,640,1280])]:
        previous=None
        for n in ns:
            label=f"g{grading:g}_n{n}"
            _,cache,meta=run_case(n,grading,rtol=1e-10,atol=1e-12,label=label)
            row={"label":label,"grading":grading,"n":n,
                 "heat_bessel_max_abs":float(np.max(np.abs(cache['temperature'][1:]-ref512[1:]))),
                 "heat_balance":meta['fields']['temperature']['independent_balance_max_abs'],
                 "water_balance":meta['fields']['moisture']['independent_balance_max_abs']}
            if previous is not None:
                for field in ['temperature','moisture']:
                    delta=np.abs(cache[field][1:]-previous[field][1:])
                    where=np.unravel_index(np.argmax(delta),delta.shape)
                    row[field+'_grid_delta_max']=float(delta[where])
                    row[field+'_grid_delta_location']={'time_s':int(where[0]+1),'radius_cm':round(float(radii[where[1]]*100),5)}
                    row[field+'_rounding_differences']=int(np.count_nonzero(np.round(cache[field][1:],4)!=np.round(previous[field][1:],4)))
            print(json.dumps(row),flush=True)
            result.append(row)
            previous=cache
            save_json(ROOT/'reports/q1_mesh_convergence.json',{'bessel_256_512_max_abs':spectral_error,'mesh_series':result})
    np.savez_compressed(ROOT/'results/cache/q1/bessel_reference.npz',time_s=times,radius_m=radii,temperature=ref512)


if __name__=='__main__':
    main()
