import os
import sys
import colorsys
import subprocess
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbk
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import matplotlib.tri as mtri
import mpl_toolkits.mplot3d
from matplotlib import rc
import math
import numpy
import scipy
import scipy.integrate
from scipy.interpolate import splrep, splev
import copy

class ColourPallete:
    @staticmethod
    def toHex(r,g,b):
        return '#%02x%02x%02x' % (r,g,b)
    @staticmethod
    def toRGB(hexval):
        value = hexval.lstrip('#')
        lv = len(value)
        return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))
    @staticmethod
    def HSV_to_RGB(hsviter):
        #in colorsys all input coords go from 0 to 1
        #keep 0 to 1 range for hsv, but use 0-255 for RGB as usual
        rgb = colorsys.hsv_to_rgb(hsviter[0],hsviter[1],hsviter[2])

        #print "HSV_to_RGB: ",hsviter, rgb
        #return a list
        out = [1,1,1]
        for i in range(3):
            out[i] = rgb[i]*255.0
        return out
    @staticmethod
    def RGB_to_HSV(rgbiter):
        #see comment above
        hsv = colorsys.rgb_to_hsv(rgbiter[0]/255.0, rgbiter[1]/255.0, rgbiter[2]/255.0)
        out = [1,1,1]
        for i in range(3):
            out[i] = hsv[i]
        return out


    @staticmethod
    def generatePallete(nhues,startcolour=(0,1,1) ):
        #startcolour should be a tuple (h,s,v)
        #vary the hue but not the saturation or value
        sep = 1.0/nhues 
        rgbpal = []
        for i in range(nhues):
            hsv = [startcolour[0]+sep*i, startcolour[1], startcolour[2] ]
            print("HSV is ",hsv)
            rgb = ColourPallete.HSV_to_RGB(hsv)
            rgbpal.append(ColourPallete.toHex(rgb[0],rgb[1], rgb[2] ) )
        return rgbpal


    @staticmethod
    def pasteliseColour(colour, alpha):
        colorrgb = colour
        if(cbk.is_string_like(colour)==True ):
            colorrgb = ColourPallete.toRGB(colour)

        colorhsv = ColourPallete.RGB_to_HSV(colorrgb)
        newrgb = ColourPallete.HSV_to_RGB([colorhsv[0],colorhsv[1]*alpha,colorhsv[2]])
        return ColourPallete.toHex(newrgb[0],newrgb[1],newrgb[2])


    @staticmethod
    def fixedPallete(nhues):
        """Return a list of $nhues colours (as hex in rgb) to use as a pallete chosen to be distinct"""
        pallete = [
            '#ff0000',
            '#00ff00',
            '#0000ff',
            '#ff9c00',
            '#a3007f',
            '#00a000',
            '#ff00ff',
            '#00ffff',
            '#d38d4e']
        if(nhues > len(pallete)):
            raise IndexError("Not enough predefined colours in pallete\n");
        return pallete[:nhues]

    @staticmethod
    def palleteMap(keylist, sort_cmp = None, reverse = False):
        """Return a dictionary mapping each key to a colour from the static list. If sort_cmp != None, a sort
        is first performed on the keylist using the sort_cmp function. cmp puts the smallest number at the top of the list!
        Use reverse=True to reverse the list post-sort"""

        lst = keylist
    
        if(sort_cmp !=None):
            lst.sort(sort_cmp)
            if(reverse==True):
                lst.reverse()

        out = {}
        pidx = 0
        for i in range(len(lst)):
            if(i>0):
                #check for multiple entries with same value, assign same colour
                if(sort_cmp != None and sort_cmp(lst[i],lst[i-1])==0):
                    print(lst[i], lst[i-1])
                    out[lst[i]] = out[lst[i-1]]
                elif(sort_cmp == None and lst[i] == lst[i-1]):
                    out[lst[i]] = out[lst[i-1]]
                else:
                    out[lst[i]] = pidx
                    pidx+=1
            else:
                out[lst[i]] = pidx
                pidx+=1

        pallete = ColourPallete.fixedPallete(pidx)        
        for i in range(len(lst)):
            out[lst[i]] = pallete[out[lst[i] ] ]

        return out

class DataSet:
        def __init__(self):
                self.x = None
                self.y = None
                self.dxm = None
                self.dxp = None
                self.dym = None
                self.dyp = None
        def createEmpty(self): 
            self.x = []
            self.y = []
            self.dxp = []
            self.dxm = []
            self.dyp = []
            self.dym = []

        def shiftX(self, shift):
            for i in range(len(self.x)):
                self.x[i] = self.x[i] + shift

        def setNaN(self, val=0.):
            for i in range(len(self.y)):
                if(math.isnan(float(self.y[i]))):
                    self.y[i] = val
            for i in range(len(self.dyp)):
                if(math.isnan(float(self.dyp[i]))):
                    self.dyp[i] = val
            for i in range(len(self.dym)):
                if(math.isnan(float(self.dym[i]))):
                    self.dym[i] = val
            
        def append(self, x,y,dx,dy):
            if(self.x == None):
                self.createEmpty()

            self.x.append(float(x))
            self.y.append(float(y))
            self.dxm.append(float(dx))
            self.dxp.append(float(dx))
            self.dym.append(float(dy))
            self.dyp.append(float(dy))

        def printit(self):
            if(self.x != None):
                for i in range(len(self.x)):
                    print("%f %f %f %f" % (self.x[i],self.y[i],self.dxp[i],self.dyp[i]))


class ErrorBand:
        def __init__(self):
                self.x = None
                self.upper = None
                self.lower = None

        def createEmpty(self): 
            self.x = []
            self.upper = []
            self.lower = []

        #Convert a DataSet into an ErrorBand. Assumes no x errors
        def importDataSet(self, dset):
            if(self.x == None):
                self.createEmpty()

            self.x = copy.deepcopy(dset.x)
            for i in range(len(dset.y)):
                self.upper.append(dset.y[i] + dset.dyp[i])
                self.lower.append(dset.y[i] - dset.dym[i])
                

class HistogramData:
        def __init__(self):
                self.y = None
                
class DataSet3D:
        def __init__(self):
                self.x = None
                self.y = None
                self.z = None

class ErrorLine:
    def __init__(self):
        self.start = None   #array of tuples of coords of line start
        self.end = None #array of tuples of coords of line end

#For data sets where the x and y coordinates are correlated but in a non-linear fashion, we can draw a curve through successive x,y points
class ErrorCurve:
    def __init__(self):
        self.curves = None   #array of arrays of tuples  [ [(x0,y0),(x1,y1)... ],   [(x0,y0),(x1,y1)... ] ... ]
        self.markers = None  #Optional locations of markers (eg at middle of curve) as array of tuples

#ecurve should be an instance of ErrorCurve
def plotErrorCurves(axes, ecurve, **kwargs):
    marker = "o"
    if 'marker' in kwargs.keys():
        marker = kwargs['marker']
        del kwargs['marker']

    if 'color' not in kwargs.keys():
        kwargs['color'] = "r"    
    if 'mec' not in kwargs.keys():
        kwargs['mec'] = 'black'
    if 'ms' not in kwargs.keys():
        kwargs["ms"] = 5 #marker size

    plot_linestyle = ''
    lines_linestyle = '-'
    if 'plot_linestyle' in kwargs.keys():
        plot_linestyle = kwargs['plot_linestyle']
        del kwargs['plot_linestyle']
    if 'lines_linestyle' in kwargs.keys():
        lines_linestyle = kwargs['lines_linestyle']
        del kwargs['lines_linestyle']


    smooth = False
    if 'smooth' in kwargs.keys():
        smooth = kwargs['smooth']
        del kwargs['smooth']
    smooth_k = 3  #Choose the order of the spline interpolation. Note: if the amount of data is <= to this number, k will be reduced to ndata-1
    if 'smooth_k' in kwargs.keys():
        smooth_k = kwargs['smooth_k']
        del kwargs['smooth_k']
    
    kwargs['linestyle'] = lines_linestyle

    lineset = []
    centerx = []
    centery = []

    for curve in ecurve.curves:
        xdata = []
        ydata = []
        for c in curve:
            xdata.append(c[0])
            ydata.append(c[1])

        if smooth == True:
            k=3
            if len(ydata) <= k:
                k = len(ydata)-1

            tck = splrep(xdata, ydata, k=k)
            xnew = numpy.linspace(xdata[0], xdata[-1])
            ynew = splev(xnew, tck)
            xdata = xnew
            ydata = ynew

        line = lines.Line2D(xdata,ydata,**kwargs)
        lineset.append(line)
        axes.add_line(line)

    kwargs['marker'] = marker
    kwargs['linestyle'] = plot_linestyle

    markers = None
    if ecurve.markers != None:
        xdata = []
        ydata = []
        for m in ecurve.markers:
            xdata.append(m[0])
            ydata.append(m[1])
        markers = axes.plot(xdata,ydata,**kwargs)
    return (lineset,markers)





#eline should be an instance of ErrorLine
def plotErrorLines(axes, eline, **kwargs):
    npt = len(eline.start)
    if(len(eline.end) != npt):
        print("plotErrorLines Error: start and end must be arrays of same length")
        exit

    marker = "o"
    if 'marker' in kwargs.keys():
        marker = kwargs['marker']
        del kwargs['marker']

    if 'color' not in kwargs.keys():
        kwargs['color'] = "r"    
    if 'mec' not in kwargs.keys():
        kwargs['mec'] = 'black'
    if 'ms' not in kwargs.keys():
        kwargs["ms"] = 5 #marker size

    plot_linestyle = ''
    lines_linestyle = '-'
    if 'plot_linestyle' in kwargs.keys():
        plot_linestyle = kwargs['plot_linestyle']
        del kwargs['plot_linestyle']
    if 'lines_linestyle' in kwargs.keys():
        lines_linestyle = kwargs['lines_linestyle']
        del kwargs['lines_linestyle']

    kwargs['linestyle'] = lines_linestyle

    lineset = []
    centerx = []
    centery = []

    for i in range(npt):
        line = lines.Line2D([eline.start[i][0],eline.end[i][0]], [eline.start[i][1],eline.end[i][1]],**kwargs)
        lineset.append(line)
        axes.add_line(line)

        #Add a marker at the center
        centerx.append(eline.start[i][0] + (eline.end[i][0]-eline.start[i][0])/2)
        centery.append(eline.start[i][1] + (eline.end[i][1]-eline.start[i][1])/2)

    kwargs['marker'] = marker
    kwargs['linestyle'] = plot_linestyle

    markers = axes.plot(centerx,centery,**kwargs)
    return (lineset,markers)

def plotDataSet(axes,dataset, **kwargs):
    if 'linestyle' not in kwargs.keys():
        kwargs['linestyle'] = ""
    if 'marker' not in kwargs.keys():
        kwargs['marker'] = "o"
    if 'ms' not in kwargs.keys():
        kwargs["ms"] = 5 #marker size

    #My extra arguments
    capwidth = 2.0
    bar_linestyle = 'solid' #ACCEPTS: ['solid' | 'dashed', 'dashdot', 'dotted' | (offset, on-off-dash-seq) ]
    color = 'r'
    hollowsymbol = False
 
    if 'capwidth' in kwargs.keys():
        capwidth = kwargs['capwidth']
        del kwargs['capwidth']
    if 'bar_linestyle' in kwargs.keys():
        bar_linestyle = kwargs['bar_linestyle']
        del kwargs['bar_linestyle']
    if 'hollowsymbol' in kwargs.keys():
        hollowsymbol = kwargs['hollowsymbol']
        del kwargs['hollowsymbol']
    if 'color' in kwargs.keys():
        color = kwargs['color']
        #del kwargs['color']

    if 'ecolor' not in kwargs.keys():
        kwargs['ecolor'] = color
    if 'mfc' not in kwargs.keys():
        kwargs['mfc'] = color
    if 'mec' not in kwargs.keys():
        kwargs['mec'] = 'black'


    if(len(dataset.x)==0):
        return

    x = dataset.x
    y = dataset.y
    dx = [dataset.dxm,dataset.dxp]
    dy = [dataset.dym,dataset.dyp]

    plotset = None
    try:
        plotset = axes.errorbar(x,y,xerr=dx,yerr=dy,**kwargs)
    except Exception as Message:
        raise(Exception, "Error plotting dataset %s" % (Message))

    #bar linestyle
    plotset[2][0].set_linestyles(bar_linestyle)
    plotset[2][1].set_linestyles(bar_linestyle)

    #capsize
    for cap in plotset[1]:
        plt.setp(cap,mew=capwidth)
    
    #Hollow symbols if required
    if(hollowsymbol == True):
        s = plotset[0]
        s.set_markerfacecolor("None")
        s.set_markeredgecolor(color)
        if 'markeredgewidth' not in kwargs.keys():
            s.set_markeredgewidth(1.25)

    return plotset



def plotErrorBand(axes, band, **kwargs):
    x = band.x
    uy = band.upper
    ly = band.lower

    usetransparency=False
    boundary_lines = False
    boundary_lines_zorder = 10

    if "usetransparency" in kwargs.keys():
        usetransparency = kwargs["usetransparency"]
        del kwargs["usetransparency"]
    if "boundary_lines" in kwargs.keys():
        boundary_lines = kwargs["boundary_lines"]
        del kwargs["boundary_lines"]
    if "boundary_lines_zorder" in kwargs.keys():
        boundary_lines_zorder = kwargs["boundary_lines_zorder"]
        del kwargs["boundary_lines_zorder"]
    if "color" not in kwargs.keys():
        kwargs["color"] = 'r'
        
    plot_band = axes.fill_between(x,ly,uy,**kwargs)

    if(boundary_lines == True):
        colour = plt.getp(plot_band,'facecolors')
        chex = ColourPallete.toHex(colour[0][0]*255,colour[0][1]*255,colour[0][2]*255)
        axes.plot(x,uy,marker='None',linestyle='--',linewidth=1.5,color=chex,zorder=boundary_lines_zorder)
        axes.plot(x,ly,marker='None',linestyle='--',linewidth=1.5,color=chex,zorder=boundary_lines_zorder)

    return plot_band


def plotHistogram(axes, data, **kwargs):
    if "nbins" in kwargs.keys() and "bins" not in kwargs.keys():
        binwidth = (max(data.y) - min(data.y))/(kwargs["nbins"])
        kwargs["bins"]=numpy.arange(min(data.y), max(data.y) + binwidth, binwidth)
        del kwargs["nbins"]
    if "color" not in kwargs.keys():
        kwargs["color"] = 'r'
        
    plot_hist = axes.hist(data.y,**kwargs)
    return plot_hist



def plotWireframe(axes3d, data, **kwargs):
    triang = mtri.Triangulation(data.x, data.y)

    if 'color' not in kwargs.keys():
        kwargs['color'] = (0,0,0,0)
    if 'edgecolor' not in kwargs.keys():
        kwargs['edgecolor'] = 'Black'
    if 'linewidth' not in kwargs.keys():
        kwargs['linewidth'] = 0.1
    
    return axes3d.plot_trisurf(triang, data.z, **kwargs)

def plotScatter(axes3d, data, **kwargs):
    if 's' not in kwargs.keys():
        kwargs['s'] = 4  #size of marker in points^2
    if 'c' not in kwargs.keys():
        kwargs['c'] = data.z  #color is z-value dependent
    if 'cmap' not in kwargs.keys():
        kwargs['cmap'] = plt.get_cmap("seismic")  #blues to reds
    if 'linewidth' not in kwargs.keys():
        kwargs['linewidth'] = 0 #edges obscure the points usually
        
    return axes3d.scatter(data.x,data.y,data.z, **kwargs)


def getHDF5valErr(filename, elem):
    if(os.path.isfile(filename) == False):
        return "ERR"
    arg_str = "-elem \"%s\"" % elem
    out = os.popen("hdf5_print %s %s" % (filename,arg_str)).read()
    m = re.search(r'([\d\.e\+\-]+)\s\+\-\s([\d\.e\+\-]+)',out)
    if m == None:
        return "ERR"
    else:
        return (float(m.group(1)),float(m.group(2)))    
