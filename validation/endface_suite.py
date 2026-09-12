"""Reproduce the conditional end-face audit; existing matching caches reused.

Stale caches fail closed. No production one-dimensional output is overwritten.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from validation.endface_run import run
from validation import endface_kernel_checks, endface_summary, endface_diagnostics
from validation import q1_nonlinear_independent, unified_error_budget


def main():
    cases=[]
    for group in ("q1","q23","q4"):
        for nr in (32,64,128):cases.append((group,nr,3,0.,1e-8,120.))
        for nr,nz in ((32,32),(64,32),(64,64),(128,64)):
            cases.append((group,nr,nz,1.,1e-8,120.))
    for group in ("q23","q4"):cases.append((group,128,128,1.,1e-8,120.))
    for group,nr,rt,dt in (("q23",64,1e-9,60.),("q4",64,1e-10,30.),("q4",128,1e-10,30.)):
        cases.extend(((group,nr,3,0.,rt,dt),(group,nr,64,1.,rt,dt)))
    assert len(cases)==29 and len(set(cases))==29
    for group,nr,nz,em,rt,dt in cases:
        _,meta=run(group,nr,nz,em,rt,dt)
        assert meta["water_balance_relative"]<1e-6
        print("verified",meta["label"],flush=True)
    q1_nonlinear_independent.main()
    endface_kernel_checks.main()
    endface_summary.main()
    endface_diagnostics.main()
    unified_error_budget.main()


if __name__=="__main__":main()
