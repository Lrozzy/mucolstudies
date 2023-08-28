#!/usr/bin/env python
# This is based on this report
# https://arxiv.org/abs/2109.00845
# https://inspirehep.net/literature/1915353
# which gives a nice overview of the different helix parameters.

import os
import math
import sys
import numpy as np
import ROOT
import matplotlib
from matplotlib import pyplot

# Set the global ROOT style to something close to the official ATLAS style.
ROOT.gROOT.SetStyle("ATLAS")

# I use these macros to set up the ATLAS root style, but unless you have access to them
# this obviously won't work. The built-in ATLAS Style (set above) is close enough.
#ROOT.gROOT.LoadMacro("AtlasLabels.C")
#ROOT.gROOT.LoadMacro("AtlasStyle.C")
#ROOT.SetAtlasStyle()

# perigee parametrization for an example track.
# I literally just made up these numbers.
d0 = 0
z0 = 0
eta = 0.3
phi = 0.3
# q/pT in 1/GeV. So this would be a 1 GeV track.
qoverpt = 1.0
# The bfield, in Tesla.
bfield = 5

# the x/y range to draw.
xmin = -1000
xmax = 1000
ymin = -1000
ymax = 1000

def getTransverseImpact(d0 = 0, phi, xref=0, yref=0):
    """ Given d0 and phi, figure out what x0 and y0 should be.
        to do this one needs to know the reference point, which is assumed to be the origin."""
    x0 = xref - d0 * np.sin(phi)
    y0 = yref + d0 * np.cos(phi)
    return x0, y0

def getCurvature(qoverpt, bfield):
    """ Convert q/pt -> curvature (1/radius of curvature). Requires B field."""
    c = .29979 * bfield * qoverpt
    return c

def getX(arclength, x0, c, phi):
    """ Given an arclength, as well as the helix parameters, returns the x coordinate at that location."""
    x = x0 + (np.sin(c*arclength/2)) / (c/2) * np.cos(phi - c*arclength/2)
    return x

def getY(arclength, y0, c, phi):
    """ Given an arclength, as well as the helix parameters, returns the y coordinate at that location."""
    y = y0 + (np.sin(c*arclength/2)) / (c/2) * np.sin(phi - c*arclength/2)
    return y

def main():
    # Convert (d0, phi, q/pt) -> (x0, y0, C).
    x0, y0 = getTransverseImpact(d0, phi)
    c = getCurvature(qoverpt, bfield)

    # run from 0 to half the arclength.
    arclengths = np.linspace(0, math.pi*(1/c))
    
    # Get a list of the xcoords and ycoords.
    xcoords = getX(arclengths, x0, c, phi) * 1000
    ycoords = getY(arclengths, y0, c, phi) * 1000
    
    # This would draw using matplotlib.
    #pyplot.plot(xcoords, ycoords)
    #pyplot.show()
    #input()
    
    # Here's a potential way to do this with ROOT.
    canvas = ROOT.TCanvas("track", "track", 800, 600)
    canvas.cd()
    
    # Create a TGraph (we could also use a 2D histogram but the TGraph draw options are a bit better).
    graph = ROOT.TGraph(len(xcoords), xcoords, ycoords)
    
    # See TGraphPainter for the meaning of these options.
    # "A" means "draw axes" (if you omit it it's equivalent to the "same" drawoption for histograms).,
    # "C" means to connect the points together with a line.
    graph.Draw("AC")
    
    # graphs are annoying to work with
    graph.GetXaxis().SetLimits(xmin, xmax)
    graph.GetYaxis().SetRangeUser(ymin, ymax)
    graph.GetXaxis().SetTitle("x [mm]")
    graph.GetYaxis().SetTitle("y [mm]")
    
    # Save the output image.
    # ROOT save PNG images by default with horrible resolution, it's usually better
    # to use a vector graphics format and then convert manually to PNG if you want that,
    # since vector graphics have "infinite" resolution/precision.
    canvas.SaveAs("track.eps")
    
if __name__ == '__main__':
    main()