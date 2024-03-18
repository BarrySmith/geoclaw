"""
Make piecewise linear topography for wave tank.
"""

from pylab import *
from clawpack.geoclaw_1d import nonuniform_grid_tools


mx = 8000  # desired number of grid cells 
print('mx = %i' % mx)

x0 = 0.
x1 = 80e3

z0_ocean = -4000.     # depth of ocean = depth at x0_slope

xzpairs = [(x0, z0_ocean),
           (x1, z0_ocean)]
           
topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

hmin = 50.  # use uniform grid in shallower water

nonuniform_grid_tools.make_celledges_cfl(x0, x1, mx, topo_fcn,
            hmin, fname='celledges.data', plot_topo=True)

