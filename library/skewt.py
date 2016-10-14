#Copyright (c) 2008 Ryan May

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.ticker import FixedLocator, AutoLocator, ScalarFormatter
import matplotlib.transforms as transforms
import matplotlib.axis as maxis
import matplotlib.artist as artist
from matplotlib.projections import register_projection

import numpy as np

import pylab

import matplotlib.pyplot as plt

#TODO:
#   *Panning and zooming are horribly broken, probably because the
#    skewed data are used as bounds. Needs to be disabled (especially panning)
#    or updated to work sensibly
#   *How do we automatically pick appropriate limits so that all relevant
#    slanted gridlines are added?
#   *How do we get the labels we want at the top?
#   *Set good aspect ratio, at least by default
#    - set_aspect(1/11.17, adjustable='datalim') would seem to do what I want
#       except that when I resize the figure, there's a lot of jumping around
#       Might be related to panning/zooming problems
#   *New functions/methods needed to add the various significant lines:
#       -Moist adiabats: ThetaE (HOW THE HECK DOES IT WORK?)
#   *Combine skewT plot with separate subplot for vertical wind plot using
#    barbs?

class SkewXTick(maxis.XTick):
    def draw(self, renderer):
        if not self.get_visible(): return
        renderer.open_group(self.__name__)

        if self.gridOn:
            self.gridline.draw(renderer)
        if self.tick1On:
            self.tick1line.draw(renderer)
        if self.tick2On:
            self.tick2line.draw(renderer)

        if self.label1On:
            self.label1.draw(renderer)
        if self.label2On:
            self.label2.draw(renderer)

        renderer.close_group(self.__name__)

    def set_clip_path(self, clippath, transform=None):
        artist.Artist.set_clip_path(self, clippath, transform)
        self.tick1line.set_clip_path(clippath, transform)
        self.tick2line.set_clip_path(clippath, transform)
        self.gridline.set_clip_path(clippath, transform)
    set_clip_path.__doc__ = artist.Artist.set_clip_path.__doc__

class SkewXAxis(maxis.XAxis):
    def _get_tick(self, major):
        return SkewXTick(self.axes, 0, '', major=major)

    def draw(self, renderer, *args, **kwargs):
        'Draw the axis lines, grid lines, tick lines and labels'
        ticklabelBoxes = []
        ticklabelBoxes2 = []

        if not self.get_visible(): return
        renderer.open_group(__name__)
        interval = self.get_view_interval()
        for tick, loc, label in self.iter_ticks():
            if tick is None: continue
            if transforms.interval_contains(interval, loc):
                tick.set_label1(label)
                tick.set_label2(label)
            tick.update_position(loc)
            tick.draw(renderer)
            if tick.label1On and tick.label1.get_visible():
                extent = tick.label1.get_window_extent(renderer)
                ticklabelBoxes.append(extent)
            if tick.label2On and tick.label2.get_visible():
                extent = tick.label2.get_window_extent(renderer)
                ticklabelBoxes2.append(extent)

        # scale up the axis label box to also find the neighbors, not
        # just the tick labels that actually overlap note we need a
        # *copy* of the axis label box because we don't wan't to scale
        # the actual bbox

        self._update_label_position(ticklabelBoxes, ticklabelBoxes2)

        self.label.draw(renderer)

        self._update_offset_text_position(ticklabelBoxes, ticklabelBoxes2)
        self.offsetText.set_text( self.major.formatter.get_offset() )
        self.offsetText.draw(renderer)

class SkewXAxes(Axes):
    # The projection must specify a name.  This will be used be the
    # user to select the projection, i.e. ``subplot(111,
    # projection='skewx')``.
    name = 'skewx'

    def _init_axis(self):
        #Taken from Axes and modified to use our modified X-axis
        "move this out of __init__ because non-separable axes don't use it"
        self.xaxis = SkewXAxis(self)
        self.yaxis = maxis.YAxis(self)
        self._update_transScale()

#    def get_axes_patch(self):
#        """
#        Override this method to define the shape that is used for the
#        background of the plot.  It should be a subclass of Patch.

#        In this case, it is a Circle (that may be warped by the axes
#        transform into an ellipse).  Any data and gridlines will be
#        clipped to this shape.
#        """
#        return Circle((0.5, 0.5), 0.5)

    def draw(self, *args):
        '''
        draw() is overridden here to allow the data transform to be updated
        before calling the Axes.draw() method.  This allows resizes to be
        properly handled without registering callbacks.  The amount of
        work done here is kept to a minimum.
        '''
        self._update_data_transform()
        Axes.draw(self, *args)

    def _update_data_transform(self):
        '''
        This separates out the creating of the data transform so that
        it alone is updated at draw time.
        '''
        # This transforms x in pixel space to be x + the offset in y from
        # the lower left corner - producing an x-axis sloped 45 degrees
        # down, or x-axis grid lines sloped 45 degrees to the right
        self.transProjection.set(transforms.Affine2D(
            np.array([[1, 1, -self.bbox.ymin], [0, 1, 0], [0, 0, 1]])))

        # Full data transform
        self.transData.set(self._transDataNonskew + self.transProjection)

    def _set_lim_and_transforms(self):
        """
        This is called once when the plot is created to set up all the
        transforms for the data, text and grids.
        """
        #Get the standard transform setup from the Axes base class
        Axes._set_lim_and_transforms(self)

        #Save the unskewed data transform for our own use when regenerating
        #the data transform. The user might want this as well
        self._transDataNonskew = self.transData

        #Create a wrapper for the data transform, so that any object that
        #grabs this transform will see an updated version when we change it
        self.transData = transforms.TransformWrapper(
            transforms.IdentityTransform())

        #Create a wrapper for the proj. transform, so that any object that
        #grabs this transform will see an updated version when we change it
        self.transProjection = transforms.TransformWrapper(
            transforms.IdentityTransform())
        self._update_data_transform()

    def get_xaxis_transform(self, which='grid'):
        """
        Get the transformation used for drawing x-axis labels, ticks
        and gridlines.  The x-direction is in data coordinates and the
        y-direction is in axis coordinates.

        We override here so that the x-axis gridlines get properly
        transformed for the skewed plot.
        """
        return self._xaxis_transform + self.transProjection

    # Disable panning until we find a way to handle the problem with
    # the projection
    def start_pan(self, x, y, button):
        pass

    def end_pan(self):
        pass

    def drag_pan(self, button, key, x, y):
        pass

# Now register the projection with matplotlib so the user can select
# it.
register_projection(SkewXAxes)

def draw_skewt(p, h, T, Td, u, v, imagename, title = None, show = False):
    '''
    If you want to plot a skewt, skewt.draw_skewt


    Inputs:
        p -> 1d numpy array, pressure in Pa/100.0 (hPa)
        h -> 1d numpy array, height in m
        t -> 1d numpy array, temp in C
        td -> 1d numpy array, dew pt temp in C
        imagename -> string, the name you want to call the image ... i.e. louise.png

    Optional Inputs:
        title -> string, obvious really
        show -> True/False, if True will open and show plot, if false, will just save to image file.

    Example usage:
        skewt.draw_skewt(p/100.0 (Grrr), h, t (C), td (C), 'louse.png', title = 'woo.png', show = False)
    
    '''
    
#    print "repr(p)"
#    print repr(p)
#    print "repr(h)"
#    print repr(h)
#    print "repr(T)"
#    print repr(T)
#    print "repr(Td)"
#    print repr(Td)
#
    plt.rcParams.update({'axes.linewidth':1,'ytick.color':'k'})

    fig = plt.figure(1, figsize=(6.5875, 6.2125))
    fig.clf()

    ax = fig.add_axes([.125,.1,.7,.8], projection='skewx')
    
    ax.grid(True)
    
    ax.semilogy(T, p, 'r')
    ax.semilogy(Td, p, 'g')
    
    xticklocs=np.arange(-80,45,10)

    # T0 = ax.get_xticks()
    T0 = xticklocs
    P0 = 1000.
    R = 287.05
    Cp = 1004.
    P = np.linspace(*ax.get_ylim()).reshape(1, -1)

    T = (T0[:,np.newaxis] + 273.15) * (P/P0)**(R/Cp) - 273.15
    linedata = [np.vstack((t[np.newaxis,:], P)).T for t in T]
    dry_adiabats = LineCollection(linedata, colors='r', linestyles='dashed',
        alpha=0.5)
    ax.add_collection(dry_adiabats)

    w = np.array([0.0004,0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024,
        0.032]).reshape(-1, 1)
    e = P * w / (0.622 + w)
    T = 243.5/(17.67/np.log(e/6.112) - 1)
    linedata = [np.vstack((t[np.newaxis,:], P)).T for t in T]
    mixing = LineCollection(linedata, colors='g', linestyles='dashed',
        alpha=0.8)
    ax.add_collection(mixing)

#    Lv = 2.4e6
#    T = T[:,0][:,np.newaxis] * (P/P0)**(R/Cp) - (Lv/Cp) * w
#    linedata = [np.vstack((t[np.newaxis,:], P)).T for t in T]
#    moist_adiabat = LineCollection(linedata, colors='b', linestyles='dashed',
#        alpha=0.8)
#    ax.add_collection(moist_adiabat)
    

    ax.set_ylabel('Pressure (hPa)')
    ax.set_xlabel('Temperature (C)')
    if isinstance(title, str):
        ax.title.set_text(title)
    
    # change plotting parameters : white box, white tick labels
    plt.rcParams.update({'axes.linewidth':0,'ytick.color':'k', 'ytick.major.size':0, 'ytick.minor.size':0})
    #plt.rcParams.update({'axes.linewidth':0,'ytick.color':'w', 'ytick.major.size':0, 'ytick.minor.size':0}) #simon changed
    
    wbax=fig.add_axes([0.85,0.1,0.125,0.8],sharey=ax)
    wbax.barbs(np.zeros(p.shape)[::2], p[::2], u[::2], v[::2])
    
    wbax.xaxis.set_ticks([],[])
    for tick in wbax.yaxis.get_major_ticks():
        tick.label1On = False

    ax.set_yticks(np.linspace(100,1000,10))
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.set_xticks(xticklocs)
    ax.set_xlim(-40,45)
    ax.set_ylim(1050,100)

    pylab.savefig(imagename)

    # if show:
    plt.draw()
    #plt.show()
