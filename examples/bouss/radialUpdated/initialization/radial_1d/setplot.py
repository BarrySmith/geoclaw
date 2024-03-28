
from pylab import *
import os, sys

try:
    from clawpack.geoclaw_1d import geoplot
except:
    print('Could not import from geoclaw_1d')


from clawpack.geoclaw_1d.nonuniform_grid_tools import make_mapc2p
from clawpack.clawutil.data import ClawData
import numpy
from importlib import reload


sys.path.insert(0,'../ic_bessel')
from radial_waves_hankel import eta_u_radial_H

import make_ic_transform
#reload(make_ic_transform)

xlimits = [0, 30e3]
true_eqn = True


h0 = make_ic_transform.h0
omega = make_ic_transform.omega
rvals = make_ic_transform.r
kvals = make_ic_transform.k
eta0hat = make_ic_transform.eta0hat
direction = make_ic_transform.direction

A = 2
g = 9.81
eta_limits = [-1.1*A, 1.1*A]
u_limits = [-1.1*A*sqrt(g/h0), 1.1*A*sqrt(g/h0)]

outdir2 = None
#outdir2 = os.path.abspath('_output_6')
#outdir2 = os.path.abspath('_output_mx5000_aneg')

def setplot(plotdata):

    plotdata.clearfigures()

    outdir1 = plotdata.outdir

    mapc2p1, mx_edge, xp_edge = make_mapc2p(os.path.join(outdir1,'celledges.txt'))

    if outdir2:
        print('Plotting from %s, comparing to %s' % (outdir1,outdir2))
    
    from clawpack.amrclaw.data import GaugeData 
    setgauges = GaugeData() 
    setgauges.read(outdir1)
    gauge_xc = {}
    for k in range(len(setgauges.gauges)):
        gauge = setgauges.gauges[k]
        gaugeno = gauge[0]
        gauge_xc[gaugeno] = gauge[1]
    

    if outdir2:
        mapc2p2, mx_edge, xp_edge = make_mapc2p(os.path.join(outdir2,'celledges.txt'))


    try:
        fname = os.path.join(plotdata.outdir, 'fgmax.txt')
        d = numpy.loadtxt(fname)
        etamax = numpy.where(d[:,1]>1e-6, d[:,3], numpy.nan)
        xmax = d[:,0]
        jmax = numpy.where(d[:,1]>0)[0].max()
        #print("run-in = %8.2f,  run-up = %8.2f" % (d[jmax,0],d[jmax,3]))
        print('Loaded hmax from ',fname)
    except:
        xmax = None
        print("Failed to load hmax from ",fname)

    xmax = None # to suppress plotting eta max

    def fixticks1(current_data):
        from pylab import ticklabel_format, grid,tight_layout
        ticklabel_format(useOffset=False)
        grid(True)
        tight_layout()
        #import pdb; pdb.set_trace()

    def fixticks(current_data):
        from pylab import ticklabel_format, plot,grid,gca
        ticklabel_format(useOffset=False)
        if xmax is not None:
            plot(xmax, etamax, 'r')
        grid(True)
        
    def add_true(current_data):
        from pylab import ticklabel_format, plot,grid,gca,exp
        ticklabel_format(useOffset=False)
        if xmax is not None:
            plot(xmax, etamax, 'r')
        grid(True)
        
        if 0:
            h0 = 4000.
            L = 400e3
            dx = L/len(current_data.x)
            x0 = 100e3
            eta0 = lambda x: exp(-((x-x0)/5e3)**2)
            t = current_data.t
            x,eta,u = eta_u(t, eta0, h0, L, dx, x0, eqn=true_eqn, make_plot=False)
            plot(x-x0,eta,'r')

        
        
    def velocity(current_data):
        from pylab import where
        q = current_data.q
        u = where(q[0,:]>1e-3, q[1,:] / q[0,:], 0.)
        return u

    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.kwargs = {'figsize':(8,6)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = eta_limits
    plotaxes.title = 'Surface displacement'
    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    if outdir2:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.surface
        plotitem.color = 'm'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p2

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = u_limits
    plotaxes.title = 'Velocity'

    def add_true(current_data):
        from pylab import ticklabel_format, plot,grid,subplot,legend
        ticklabel_format(useOffset=False)
        if xmax is not None:
            plot(xmax, etamax, 'r')
        grid(True)
        
        t = current_data.t
        
        #import pdb; pdb.set_trace()
        print('Computing transforms ...')
        eta,u = eta_u_radial_H(t,rvals,kvals,eta0hat, omega,h0,direction)
        #label = 'linearized omega'
        label = 'linearized omega (Airy)'
        subplot(211)
        plot(rvals,eta,'r',label=label)
        legend()
        subplot(212)
        plot(rvals,u,'r',label=label)
        legend()
        

    if true_eqn:
        plotaxes.afteraxes = add_true
    else:
        plotaxes.afteraxes = fixticks1

    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.velocity
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    if outdir2:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.velocity
        plotitem.color = 'm'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p2



    #----------

    plotfigure = plotdata.new_plotfigure(name='compare', figno=2)
    plotfigure.kwargs = {'figsize':(10,5)}
    plotfigure.show = (outdir2 is not None)
    

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [50e3,65e3]
    plotaxes.ylimits = [-20,50]
    plotaxes.title = 'Comparison'

    plotaxes.afteraxes = fixticks1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.outdir = outdir2
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'm'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    
    def fixgauge(current_data):
        from pylab import grid, title
        grid(True)
        gaugeno = current_data.gaugeno
        xc = gauge_xc[gaugeno]
        xp = mapc2p1(xc)
        print('+++ xc,xp:', xc,xp)
        title('Surface elevation at Gauge %i, x = %.0f m' \
              % (gaugeno, xp))

    plotaxes.afteraxes = fixgauge
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-100,150]
    plotaxes.title = 'Surface elevation eta'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 2  # eta
    plotitem.plotstyle = 'b-'


    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'      # Which frames to output
    plotdata.print_fignos = 'all'      # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output
    plotdata.parallel = True

    return plotdata

