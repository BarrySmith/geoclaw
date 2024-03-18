
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 


import numpy
from numpy import sqrt,exp,linspace
import glob, os, sys
from importlib import reload
from numpy import loadtxt
from matplotlib import image

#from clawpack.geoclaw_1d.nonuniform_grid_tools import make_mapc2p
from nonuniform_grid_tools import make_mapc2p  # local version

import matplotlib
matplotlib.rcParams['animation.embed_limit'] = 2**128

try:
    BoussDev = os.environ['BoussDev']
except:
    BoussDev = '/mnt/home/bsmith/clawpack/geoclaw/examples/bouss/'
    print("*** Need to define environment variable BoussDev to path")

new_python_dir = os.path.join(BoussDev, 'new_python')
sys.path.insert(0, new_python_dir)
import gridtools
sys.path = sys.path[1:]  # remove from path

sys.path.insert(0,BoussDev+'/initialization/ic_bessel')
#from make_ic import eta_u_radial, H1transform
from radial_waves_hankel import eta_u_radial_H, Htransform

try:
    reload(make_ic_bessel)  # in case it was changed
except:
    import make_ic_bessel

h0 = 4000.
g = 9.81
etamax = 1.

h0 = make_ic_bessel.h0
g = 9.81
omega = make_ic_bessel.omega
rvals = make_ic_bessel.r
kvals = make_ic_bessel.k
direction = make_ic_bessel.direction
print('make_ic_bessel.direction = ',direction)

r_eta_u = loadtxt('starting.data',skiprows=2)
r = r_eta_u[:,0]
eta0 = r_eta_u[:,1]
u0 = r_eta_u[:,2]

#ketahat = loadtxt('etahat.data',dtype=numpy.complex128)
#kvals = real(ketahat[:,0])
#etahat1 = ketahat[:,1]

etahat1 = Htransform(r,eta0,kvals)


def add_true(current_data):
    from pylab import ticklabel_format, plot,grid,subplot
    ticklabel_format(useOffset=False)
    grid(True)
    
    t = current_data.t
    
    #import pdb; pdb.set_trace()
    print('Computing transforms ...')
    eta,u = eta_u_radial_H(t,r,kvals,etahat1, omega,h0,direction=direction)
    subplot(211)
    plot(r,eta,'r')
    subplot(212)
    plot(r,u,'r')


#outdir_1d = None
outdir_1d = os.path.abspath('1d_radial/_output')

#outdir_1d = '1d_radial_multiphys/_output'
if outdir_1d:
    print('Comparing to 1d solution in ', outdir_1d)
    celledges_file = os.path.join(outdir_1d,'celledges.data')
    mapc2p1, mx_edge, xp_edge = make_mapc2p(celledges_file)
else:
    mapc2p1 = None



# --------------------------
def setplot(plotdata=None):
# --------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items dat
    plotdata.format = "ascii"
    #plotdata.format = "binary"

    try:
        tsudata = open(plotdata.outdir+'/geoclaw.data').readlines()
        for line in tsudata:
            if 'sea_level' in line:
                sea_level = float(line.split()[0])
                print("sea_level = ",sea_level)
    except:
        print("Could not read sea_level, setting to 0.")
        sea_level = 0.



    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             format_string='k.', add_labels=True, fontsize=6)
             #gaugenos=range(1000,1067), format_string='k.', add_labels=True,

    def timeformat(t):
        from numpy import mod
        hours = int(t/3600.)
        tmin = mod(t,3600.)
        min = int(tmin/60.)
        sec = int(mod(tmin,60.))
        timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
        return timestr
        
    def title_hours(current_data):
        from pylab import title
        t = current_data.t
        timestr = timeformat(t)
        title('Eta at t = %.0f seconds' % t, fontsize=16)

    def aframe(current_data):
        from pylab import figure, savefig

    

    def surface_or_depth(current_data):
        """
        Modified from geoplot version to use eta = q[-1,:,:], which
        should work for either num_eqn = 3 or 5.
        
        Return a masked array containing the surface elevation where the topo is
        below sea level or the water depth where the topo is above sea level.
        Mask out dry cells.  Assumes sea level is at topo=0.
        Surface is eta = h+topo, assumed to be output as 4th column of fort.q
        files.
        """
        import numpy

        #drytol = getattr(current_data.user, 'drytol', drytol_default)
        drytol = 1e-3
        q = current_data.q
        h = q[0,:,:]
        eta = q[-1,:,:]
        topo = eta - h

        # With this version, the land is transparent.
        surface_or_depth = numpy.ma.masked_where(h <= drytol,
                                                 numpy.where(topo<0, eta, h))

        try:
            # Use mask covering coarse regions if it's set:
            m = current_data.mask_coarse
            surface_or_depth = numpy.ma.masked_where(m, surface_or_depth)
        except:
            pass

        return surface_or_depth

    def topo(current_data):
       """
       Return topography = eta - h.
       Surface eta is assumed to be output as last column of fort.q files.
       Modified from geoplot version.
       """
       q = current_data.q
       h = q[0,:,:]
       eta = q[-1,:,:]
       topo = eta - h
       return topo

    #-----------------------------------------
    # Figure for big area
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Pacific', figno=0)
    #plotfigure.kwargs = {'figsize': (9,10)}
    plotfigure.show = False
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.cmd = 'subplot(121)'
    plotaxes.title = 'Pacific'
    plotaxes.scaled = False
    # modified for debug:
    plotaxes.xlimits = [-129,-122]
    plotaxes.ylimits = [43,50]

    def aa(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi, savefig
        title_hours(current_data)
        ticklabel_format(useOffset=False)
        xticks(rotation=20)
        a = gca()
        #a.set_aspect(1./cos(46.86*pi/180.))
        #addgauges(current_data)
    plotaxes.afteraxes = aa

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = surface_or_depth  # local version
    my_cmap = colormaps.make_colormap({-1.0: [0.0,0.0,1.0], \
                                     -0.5: [0.5,0.5,1.0], \
                                      0.0: [1.0,1.0,1.0], \
                                      0.5: [1.0,0.5,0.5], \
                                      1.0: [1.0,0.0,0.0]})
    plotitem.imshow_cmap = my_cmap
    #plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -0.3
    plotitem.imshow_cmax = 0.3
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = topo
    from numpy import arange, linspace
    plotitem.contour_levels = linspace(-6000,0,7)
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,1,0]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(0., 11., 1.)
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,0,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for surface with transect too
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=20)
    plotfigure.kwargs = {'figsize': (8,8)}
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.xlimits = [-50e3,50e3]
    plotaxes.ylimits = [-50e3,50e3]
    plotaxes.afteraxes = aa
    plotaxes.axescmd = 'axes([.15,.5,.7,.45])'

    # Water
    #plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = surface_or_depth  # local version
    plotitem.imshow_cmap = my_cmap
    #plotitem.imshow_cmap = geoplot.tsunami_colormap
    #plotitem.imshow_cmin = -0.3 * etamax
    #plotitem.imshow_cmax = 0.3 * etamax
    plotitem.imshow_cmin = -0.1 * etamax
    plotitem.imshow_cmax = 0.1 * etamax
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0,1,1,1,1,1]
    plotitem.amr_patchedges_color = ['k','r','g','c','b','m']


    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = topo
    plotitem.contour_levels = [-2900, -2000, -1000, -201, 0]
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'-','linewidths':1}
    plotitem.amr_contour_show = [0,0,0,0,1,0]
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0,0,0,0,1,0]



    #-----------------------------------------
    # Figure for cross section compared to 1d_radial
    #-----------------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'Transect of eta'

    def plot_xsec(current_data):
        from pylab import plot,legend,xlabel,sqrt,grid,xlim,ylim,title
        from numpy import cos,pi,linspace,zeros,ones
        from clawpack.pyclaw import Solution
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
        eps = 1e-3
        # next 3 lines constant in y over all x
        xout = linspace(-50e3,50e3,2001)
        ytrans = 1.
        yout = ytrans*ones(xout.shape) + eps
        # next 3 lines constant in x over all y
        #yout = linspace(-50e3,50e3,2001)
        #xtrans = 1.
        #xout = xtrans*ones(yout.shape) + eps
        # next 3 lines along diagonal x = y. WRONG
        #yout = linspace(-50e3,50e3,2001) 
        #xout = yout 
        #rout = sqrt(xout**2 + yout**2 )
        
        if 0:
            eta_out = gridtools.grid_output_2d(framesoln, -1, xout, yout)
            #h_out = gridtools.grid_output_2d(framesoln, 0, xout, yout)
            #B_out = eta_out - h_out
            plot(xout, eta_out, 'b', label='along y=%.3f' % ytrans)

            eta_out = gridtools.grid_output_2d(framesoln, -1, xout, xout)
            plot(xout*sqrt(2), eta_out, 'g', label='along x=y')
        else:
        #elif 0:
            #method = 'nearest' # the default, which gives pw constant
            method = 'linear'  # pw linear interpolation
            eta_out1 = gridtools.grid_output_2d(framesoln, -1, xout, yout,
                                               levels=[1],method=method)
            eta_out2 = gridtools.grid_output_2d(framesoln, -1, xout, yout,
                                               levels=[2],method=method)
            eta_out3 = gridtools.grid_output_2d(framesoln, -1, xout, yout,
                                               levels=[3],method=method)
            eta_out4 = gridtools.grid_output_2d(framesoln, -1, xout, yout,
                                               levels=[4],method=method)
            eta_out5 = gridtools.grid_output_2d(framesoln, -1, xout, yout,
                                               levels=[5],method=method)
            eta_out = gridtools.grid_output_2d(framesoln, -1, xout, yout,
                                               method=method)
            #plot(xout, eta_out3, 'g', label='2D level 3')
            #plot(xout, eta_out4, 'c', label='2D level 4')
            #plot(xout, eta_out5, 'b', label='2D level 5')

            # for horizontal line  plot
            plot(xout, eta_out, 'm', label='2D transect')
            # for vertical line  plot
            #plot(yout, eta_out, 'm', label='2D transect')
            # for diagonal plot, against radial distance
            #plot(rout, eta_out, 'm', label='2D transect')
            
            framesoln = Solution(frameno, path=outdir_1d, file_format='ascii')
            xc_1d = framesoln.patch.grid.c_centers[0]
            xp_1d = mapc2p1(xc_1d)
            s = framesoln.states[0]  # only one grid patch
            q = s.get_q_global()
            eta_1d = q[-1,:]
            plot(xp_1d, eta_1d, 'k--', label='1D radial')
            plot(-xp_1d, eta_1d, 'k--')

        xlabel('radial distance',fontsize=12)
        grid(True)

        t = current_data.t
        if 0:
            print('Computing transforms ...')
            eta,u = eta_u_radial_H(t,r,kvals,etahat1, omega,h0,direction=direction)
            #plot(r, eta0, 'k', label='eta0')
            plot(-r, eta, 'r')
            xlim(-75e3,75e3)
        if 0:
            xmax = 10e3 + sqrt(g*h0)*t
            xlim(-xmax,xmax)
        #ylim(-40,80)
        legend(loc='upper left', fontsize=10)
        t = current_data.t
        title('Transect of eta at t = %.0f seconds' % t, fontsize=16)
        xlim(-50e3,50e3)

    plotaxes.afteraxes = plot_xsec

    def eta(current_data):
        q = current_data.q
        eta = q[-1,:]
        return eta

    #if outdir_1d:
    if 0:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir_1d
        plotitem.plot_var = eta
        plotitem.plotstyle = 'k-'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p1


    #-----------------------------------------
    # Velocity
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='Velocity_H', figno=21)
    plotfigure.kwargs = {'figsize': (8,8)}
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.xlimits = [-50e3,50e3]
    plotaxes.ylimits = [-50e3,50e3]
    plotaxes.afteraxes = aa
    plotaxes.axescmd = 'axes([.15,.5,.7,.45])'

    def urad(current_data):
        q = current_data.q
        u = q[1,:,:]/q[0,:,:]
        v = q[2,:,:]/q[0,:,:]
        x = current_data.x
        y = current_data.y
        urad = (u*x + v*y) / sqrt(x**2 + y**2 + 1e-4)
        return urad

    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = urad
    plotitem.imshow_cmap = my_cmap
    #plotitem.imshow_cmap = geoplot.tsunami_colormap
    #plotitem.imshow_cmin = -0.01 * etamax
    #plotitem.imshow_cmax = 0.01 * etamax
    plotitem.imshow_cmin = -0.002 * etamax
    plotitem.imshow_cmax = 0.002 * etamax
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,1,1,1,1,1]
    plotitem.amr_patchedges_color = ['k','r','g','c','b','m']


    #-----------------------------------------
    # Figure for cross section compared to 1d_radial
    #-----------------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'Transect of u'

    def plot_xsec(current_data):
        from pylab import plot,legend,xlabel,sqrt,grid,xlim,ylim
        from numpy import cos,pi,linspace,zeros,ones,where,real
        from clawpack.pyclaw import Solution
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
        eps = 1e-3

        # next 3 lines for horizontal transect
        xout = linspace(-75e3,75e3,2001)
        #ytrans = 1.
        #yout = ytrans*ones(xout.shape) + eps
        ytrans = 0. + eps
        yout = ytrans*ones(xout.shape) 

        hu_out = gridtools.grid_output_2d(framesoln, 1, xout, yout)
        hv_out = gridtools.grid_output_2d(framesoln, 2, xout, yout)
        h_out = gridtools.grid_output_2d(framesoln, 0, xout, yout)
        u_out = hu_out / h_out
        u_out = where(xout < 0, -u_out, u_out)
        v_out = hv_out / h_out
        v_out = where(xout < 0, -v_out, v_out)

        
        plot(xout, u_out, 'b', label='u on y=%.3f' % ytrans)
        plot(xout, v_out, 'm', label='v on y=%.3f' % ytrans)
        #plot(xout, u_out, 'b', label='along x=%.3f' % xtrans)

        #eta_out = gridtools.grid_output_2d(framesoln, -1, xout, xout)
        #plot(xout*sqrt(2), eta_out, 'g', label='along x=y')

        xlabel('radial distance',fontsize=12)
        grid(True)

        t = current_data.t
        print('Computing transforms ...')
        eta,u = eta_u_radial_H(t,r,kvals,etahat1, omega,h0,direction=direction)
        #plot(r, u0, 'k', label='u0')
        plot(r, real(u), 'r', label='true')
        plot(-r, u, 'r')
        plot(-r, eta*sqrt(g/h0), 'c', label='evector')
        xlim(-75e3,75e3)
        if 0:
            xmax = 10e3 + sqrt(g*h0)*t
            xlim(-xmax,xmax)
        #ylim(-40,80)
        legend(loc='lower left')

    plotaxes.afteraxes = plot_xsec

    def u(current_data):
        q = current_data.q
        u = q[1,:] / q[0,:]
        return u

    #if outdir_1d:
    if 0:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir_1d
        plotitem.plot_var = u
        plotitem.plotstyle = 'k-'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p1

 
    #-----------------------------------------
    # Velocity (with vertical transect)
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='Velocity_V', figno=22)
    plotfigure.kwargs = {'figsize': (8,8)}
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.xlimits = [-50e3,50e3]
    plotaxes.ylimits = [-50e3,50e3]
    plotaxes.afteraxes = aa
    plotaxes.axescmd = 'axes([.15,.5,.7,.45])'

    def urad(current_data):
        q = current_data.q
        u = q[1,:,:]/q[0,:,:]
        v = q[2,:,:]/q[0,:,:]
        x = current_data.x
        y = current_data.y
        urad = (u*x + v*y) / sqrt(x**2 + y**2 + 1e-4)
        return urad

    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = urad
    plotitem.imshow_cmap = my_cmap
    #plotitem.imshow_cmap = geoplot.tsunami_colormap
    #plotitem.imshow_cmin = -0.03 * etamax
    #plotitem.imshow_cmax = 0.03 * etamax
    plotitem.imshow_cmin = -0.002 * etamax
    plotitem.imshow_cmax = 0.002 * etamax
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,1,1,1,1,1]
    plotitem.amr_patchedges_color = ['k','r','g','c','b','m']


    #-----------------------------------------
    # Figure for cross section compared to 1d_radial
    #-----------------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'Transect of u'

    def plot_xsec(current_data):
        from pylab import plot,legend,xlabel,sqrt,grid,xlim,ylim
        from numpy import cos,pi,linspace,zeros,ones,where,real
        from clawpack.pyclaw import Solution
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
        eps = 1e-3

        # next 3 lines for vertical transect
        yout = linspace(-75e3,75e3,2001)
        #xtrans = 1.
        xtrans = 0. + eps
        xout = xtrans*ones(yout.shape) 

        hu_out = gridtools.grid_output_2d(framesoln, 1, xout, yout)
        hv_out = gridtools.grid_output_2d(framesoln, 2, xout, yout)
        h_out = gridtools.grid_output_2d(framesoln, 0, xout, yout)
        u_out = hu_out / h_out
        u_out = where(yout < 0, -u_out, u_out)
        v_out = hv_out / h_out
        v_out = where(yout < 0, -v_out, v_out)

        
        plot(yout, u_out, 'b', label='u on x=%.3f' % xtrans)
        plot(yout, v_out, 'm', label='v on x=%.3f' % xtrans)

        #eta_out = gridtools.grid_output_2d(framesoln, -1, xout, xout)
        #plot(xout*sqrt(2), eta_out, 'g', label='along x=y')

        xlabel('radial distance')
        grid(True)

        t = current_data.t
        #print('Computing transforms ...')
        #eta,u = eta_u_radial_H(t,r,kvals,etahat1, omega,h0,direction=direction)
        ##plot(r, u0, 'k', label='u0')
        #plot(r, real(u), 'r', label='true')
        #plot(-r, u, 'r')
        #plot(-r, eta*sqrt(g/h0), 'c', label='evector')
        #xlim(-75e3,75e3)
        if 0:
            xmax = 10e3 + sqrt(g*h0)*t
            xlim(-xmax,xmax)
        #ylim(-40,80)
        legend(loc='lower left')

    plotaxes.afteraxes = plot_xsec

    def v(current_data):
        q = current_data.q
        v = q[2,:] / q[0,:]
        return u

    #if outdir_1d:
    if 0:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir_1d
        plotitem.plot_var = v
        plotitem.plotstyle = 'k-'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p1

 

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------

    def fix_gauge(current_data):
        from pylab import plot, legend, xticks, floor, yticks,\
            xlabel,savefig,xlim,where,nan,ones,grid
        t = current_data.t
        gaugeno = current_data.gaugeno
        q = current_data.q

        h = q[0,:]
        h0 = q[0,0] * ones(h.shape)
        level = current_data.gaugesoln.level
        #B = current_data.gaugesoln.aux[0,:]
        dh_refine = 0.
        for j in range(1,len(h0)):
            if level[j] != level[j-1]:
                dh_refine = dh_refine + h[j] - h[j-1]  #B[j-1]-B[j] 
            h0[j] = h0[j] + dh_refine

        ddepth = q[0,:] - h0[:]
        #plot(t, ddepth, 'b-')

        n = int(floor(t.max()/1800.) + 2)
        xticks([1800*i for i in range(n)],[str(i/2.) for i in range(n)],\
          fontsize=15)
        yticks(fontsize=15)
        xlabel("Hours")
        grid(True)
        #save_gauge(current_data)

    #-----------------------------------------
    # Momentm Corrections
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='Psi_X', figno=23)
    plotfigure.kwargs = {'figsize': (8,8)}
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.xlimits = [-50e3,50e3]
    plotaxes.ylimits = [-50e3,50e3]
    plotaxes.afteraxes = aa
    plotaxes.axescmd = 'axes([.15,.5,.7,.45])'

    def uc(current_data):
        q = current_data.q
        uc = q[3,:,:]
        vc = q[4,:,:]
        x = current_data.x
        y = current_data.y
        return uc

    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = uc
    plotitem.imshow_cmap = my_cmap
    #plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -0.0001 
    plotitem.imshow_cmax = 0.0001 
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,1,1,1,1,1]
    plotitem.amr_patchedges_color = ['k','r','g','c','b','m']


    #-----------------------------------------
    # Figure for cross section compared to 1d_radial
    #-----------------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'Horizontal Transects'

    def plot_xsec(current_data):
        from pylab import plot,legend,xlabel,sqrt,grid,xlim,ylim
        from numpy import cos,pi,linspace,zeros,ones,where,real
        from clawpack.pyclaw import Solution
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
        eps = 1e-3

        # next 3 lines for horizontal transect
        xout = linspace(-75e3,75e3,2001)
        #ytrans = 1.
        #yout = ytrans*ones(xout.shape) + eps
        ytrans = 0. + eps
        yout = ytrans*ones(xout.shape) 

        hu_out = gridtools.grid_output_2d(framesoln, 3, xout, yout)
        hv_out = gridtools.grid_output_2d(framesoln, 4, xout, yout)
        uc_out = hu_out
        #uc_out = where(xout < 0, -uc_out, uc_out)
        vc_out = hv_out 
        #vc_out = where(xout < 0, -vc_out, vc_out)

        
        plot(xout, uc_out, 'b', label='u on y=%.3f' % ytrans)
        plot(xout, vc_out, 'm', label='v on y=%.3f' % ytrans)
        #plot(xout, u_out, 'b', label='along x=%.3f' % xtrans)

        #eta_out = gridtools.grid_output_2d(framesoln, -1, xout, xout)
        #plot(xout*sqrt(2), eta_out, 'g', label='along x=y')

        xlabel('radial distance')
        grid(True)

        #t = current_data.t
        #print('Computing transforms ...')
        #eta,u = eta_u_radial_H(t,r,kvals,etahat1, omega,h0,direction=direction)
        ##plot(r, u0, 'k', label='u0')
        #plot(r, real(u), 'r', label='true')
        #plot(-r, u, 'r')
        #plot(-r, eta*sqrt(g/h0), 'c', label='evector')
        #xlim(-75e3,75e3)
        if 0:
            xmax = 10e3 + sqrt(g*h0)*t
            xlim(-xmax,xmax)
        #ylim(-40,80)
        legend(loc='lower left')

    plotaxes.afteraxes = plot_xsec

    def u(current_data):
        q = current_data.q
        u = q[1,:] / q[0,:]
        return u

    #if outdir_1d:
    if 0:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir_1d
        plotitem.plot_var = u
        plotitem.plotstyle = 'k-'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p1

    #-----------------------------------------
    # Momentm Corrections
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='Psi_Y', figno=24)
    plotfigure.kwargs = {'figsize': (8,8)}
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.xlimits = [-50e3,50e3]
    plotaxes.ylimits = [-50e3,50e3]
    plotaxes.afteraxes = aa
    plotaxes.axescmd = 'axes([.15,.5,.7,.45])'

    def vc(current_data):
        q = current_data.q
        uc = q[3,:,:]
        vc = q[4,:,:]
        x = current_data.x
        y = current_data.y
        return vc

    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = vc
    plotitem.imshow_cmap = my_cmap
    #plotitem.imshow_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmin = -0.0001 
    plotitem.imshow_cmax = 0.0001 
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,1,1,1,1,1]
    plotitem.amr_patchedges_color = ['k','r','g','c','b','m']


    #-----------------------------------------
    # Figure for cross section compared to 1d_radial
    #-----------------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'Vertical Transects'

    def plot_xsec(current_data):
        from pylab import plot,legend,xlabel,sqrt,grid,xlim,ylim
        from numpy import cos,pi,linspace,zeros,ones,where,real
        from clawpack.pyclaw import Solution
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
        eps = 1e-3

        # next 3 lines for horizontal transect
        yout = linspace(-75e3,75e3,2001)
        xtrans = 0. + eps
        xout = xtrans*ones(yout.shape) 

        hu_out = gridtools.grid_output_2d(framesoln, 3, xout, yout)
        hv_out = gridtools.grid_output_2d(framesoln, 4, xout, yout)
        uc_out = hu_out
        #uc_out = where(xout < 0, -uc_out, uc_out)
        vc_out = hv_out 
        #vc_out = where(xout < 0, -vc_out, vc_out)

        
        plot(yout, uc_out, 'b', label='u on x=%.3f' % xtrans)
        plot(yout, vc_out, 'm', label='v on x=%.3f' % xtrans)

        #eta_out = gridtools.grid_output_2d(framesoln, -1, xout, xout)
        #plot(xout*sqrt(2), eta_out, 'g', label='along x=y')

        xlabel('radial distance')
        grid(True)

        #t = current_data.t
        #print('Computing transforms ...')
        #eta,u = eta_u_radial_H(t,r,kvals,etahat1, omega,h0,direction=direction)
        ##plot(r, u0, 'k', label='u0')
        #plot(r, real(u), 'r', label='true')
        #plot(-r, u, 'r')
        #plot(-r, eta*sqrt(g/h0), 'c', label='evector')
        #xlim(-75e3,75e3)
        if 0:
            xmax = 10e3 + sqrt(g*h0)*t
            xlim(-xmax,xmax)
        #ylim(-40,80)
        legend(loc='lower left')

    plotaxes.afteraxes = plot_xsec

    def u(current_data):
        q = current_data.q
        u = q[1,:] / q[0,:]
        return u

    #if outdir_1d:
    if 0:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir_1d
        plotitem.plot_var = u
        plotitem.plotstyle = 'k-'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p1

 

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------

    def fix_gauge(current_data):
        from pylab import plot, legend, xticks, floor, yticks,\
            xlabel,savefig,xlim,where,nan,ones,grid
        t = current_data.t
        gaugeno = current_data.gaugeno
        q = current_data.q

        h = q[0,:]
        h0 = q[0,0] * ones(h.shape)
        level = current_data.gaugesoln.level
        #B = current_data.gaugesoln.aux[0,:]
        dh_refine = 0.
        for j in range(1,len(h0)):
            if level[j] != level[j-1]:
                dh_refine = dh_refine + h[j] - h[j-1]  #B[j-1]-B[j] 
            h0[j] = h0[j] + dh_refine

        ddepth = q[0,:] - h0[:]
        #plot(t, ddepth, 'b-')

        n = int(floor(t.max()/1800.) + 2)
        xticks([1800*i for i in range(n)],[str(i/2.) for i in range(n)],\
          fontsize=15)
        yticks(fontsize=15)
        xlabel("Hours")
        grid(True)
        #save_gauge(current_data)


    plotfigure = plotdata.new_plotfigure(name='gauge eta,ddepth', figno=301, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False
    plotfigure.kwargs = {'figsize':(12,5)}
    plotfigure.facecolor = 'w'


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.5*3600, 1.5*3600]
    plotaxes.ylimits = [-12,12]
    plotaxes.title = 'Surface elevation '

    #def ddepth(current_data):
    #    q = current_data.q
    #    return q[0,:] - q[0,0]
        
    # Plot ddepth as blue curve:
    #plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.plot_var = ddepth
    #plotitem.plotstyle = 'b-'

    # plot eta as red curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'r-'
    plotaxes.afteraxes = fix_gauge


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    #plotdata.print_framenos = range(24,27)           # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    #plotdata.print_fignos = [20]             # list of figures to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True

    return plotdata

    
