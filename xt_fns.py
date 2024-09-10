# X-ray transits with batman
# Copyright (C) 2024 George King
#
# Original batman code copyright (C) 2015 Laura Kreidberg
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from scipy import integrate
import astropy.units as u
from scipy import stats
from scipy.special import kn
from scipy.integrate import quad
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from astropy.constants import G

pi2 = 2*np.pi

def integrand(z, x, tau):
    """Function defining the integrand for the integration"""
    return np.exp(  -tau * np.hypot(z, x)  )

@np.vectorize
def integration(y, Rx, h):
    return quad(integrand, np.sqrt(Rx*Rx - y*y), np.inf, args=(y, 1./h) )[0]

def get_intVals(h, Nscale=6):
    """Function to generate a new model for the intensity profile of the coronal emission"""
    #define Rx (disc edge in x coordinates) from h and Nscale
    Rx = 1. - ( Nscale * h )
    #define number of steps for the integration
    Nsteps = 10001
    #set up x array for 1001 steps from x=0 to x=1
    x = np.linspace(0, 1, Nsteps)
    #two sections of the x arrays
    xIn, xOut =  x[x<=Rx], x[x>Rx]
    #calculate the integral for the current step inside photospheric disc
    wIdisc = integration(xIn, Rx, h)
    #secondly for x > Rx
    noDisc = 2*xOut*kn(1, xOut/h)
    #combine the two
    vals = np.append(wIdisc, noDisc)

    ##add in extra point at discontinuity into outside of disc arrays##
    #extra value in x array
    xOutX = np.insert(xOut, 0, Rx)
    #insert 2 times last value in wIdisc array - value for outside of disc at discontinuity
    noDiscX = np.insert( noDisc, 0, (wIdisc[-1] * 2) ) 
    #get new vals array where the boundary discontinuity value is replaced with midpoint of the two values for that point (given by the two calculations one can do at that point)
    vals[ len(xIn)-1 ] = (noDiscX[0] + wIdisc[-1]) / 2.

    #adjust values by factor to account for 2D emission greater circumference further out
    wIdisc_2pr = wIdisc * 2 * np.pi * xIn
    noDiscX_2pr = noDiscX * 2 * np.pi * xOutX

    ##Find the normalisation constant by integrating over each section and summing##
    #within disc#
    #get the difference between each x array value
    wIdx = np.diff(xIn)
    #area under curve using sum of trapezium areas
    areaIn = np.sum( 0.5 * ( wIdisc_2pr[1:] + wIdisc_2pr[:-1] ) * wIdx)
    #outside disc#
    #diff between each x value
    oTdx = np.diff(xOutX)
    #area under curve using sum of trapezium areas
    areaOut = np.sum( 0.5 * ( noDiscX_2pr[1:] + noDiscX_2pr[:-1] ) * oTdx)
    #total area under curve
    totArea = areaIn + areaOut
    
    #get the normalised values
    norm = vals / totArea
    #return the normalised values
    return norm
    
def hfromH(H, Nscale=6):
    """Converts the emission scale height in stellar radii to its value in batman coordinates"""
    #calculate h in batman coordinates
    h = H / ( 1 + (Nscale * H) )
    return h

def Hfromh(h, Nscale=6):
    """Converts the emission scale height in batman coordinates to stellar radii"""
    #calculate H in stellar radii
    H = h / ( 1 - (Nscale * h) )
    return H

def HfromT(T, Rs, Ms, mu=1.27, kT=True):
    """Function to obtain the emission scale height in stellar radii, given a temperature.
    By default, the function assumes the value given is kT in units of keV, but Kelvin can be used if kT is set to False. The mass and radius of the star should be supplied in Solar units."""
    #calc surface gravity
    g = G.value * (Ms*u.Msun.to(u.kg)) / (Rs*u.Rsun.to(u.m))**2
    #set hydrogen atom mass
    mH = 1.673e-27
    #code for if temp is in kT (units of keV)
    if kT:
        He = T*u.keV.to(u.J) / (mu * mH * g)
    else: He = (1.380649e-23 * T) / (mu * mH * g)
    #convert to stellar radii
    HeRs = He / (Rs*u.Rsun.to(u.m))
    #return the emission scale height
    return HeRs

def calcRx(h, Nscale=6):
    """Function to calcualte the photospheric edge of the star in batman coordinates"""
    return 1 - (h * Nscale)

    

