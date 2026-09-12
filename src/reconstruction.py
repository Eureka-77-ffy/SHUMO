"""Axis/surface recovery and physical-radius sampling of cell-centred states."""
import numpy as np
from scipy.interpolate import PchipInterpolator


def reconstruct(model, times, cell_values, radii, uniform_initial=None):
    m = model.mesh
    radii = np.asarray(radii)
    if np.any(radii < 0) or np.any(radii > m.radius):
        raise ValueError("Sample outside physical material")
    r0, r1 = m.centres[:2]
    axis = (r1*r1*cell_values[:,0]-r0*r0*cell_values[:,1])/(r1*r1-r0*r0)
    surface = model.boundary(times,cell_values[:,-1])
    points = np.r_[0.,m.centres,m.radius]
    values = np.column_stack((axis,cell_values,surface))
    sampled = PchipInterpolator(points,values,axis=1,extrapolate=False)(radii)
    near_axis = radii < r0
    # Even quadratic on the central interval enforces derivative zero at the axis.
    if np.any(near_axis):
        curvature = (cell_values[:,1]-cell_values[:,0])/(r1*r1-r0*r0)
        sampled[:,near_axis] = axis[:,None]+curvature[:,None]*radii[near_axis]**2
    if uniform_initial is not None:
        sampled[np.asarray(times)==0]=uniform_initial
    return sampled
