#! /usr/bin/env python
#######################################################################
##  Michael Gordon
##  June 2014
##  Computer Club 2014 / HW1
##    Spectral feature analyzer
##
##  usage: gordon.py [-h] file
##
##  Calculate line ratios of input file.
##
##  positional arguments:
##   file        Input file for analysis
##
##  optional arguments:
##    -h, --help  show this help message and exit
#######################################################################
import numpy as np
import os
import pyfits
import argparse
from scipy.optimize import curve_fit
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from scipy.special import wofz

# Fitting functions
def gaussian(x, *p):
    A, mu, sigma, offset = p
    return A*np.exp(-(x-mu)**2/(2.0*sigma**2))+offset

def voigt(x, *p):
    A, mu, sigma, gamma = p
    z = ((x-mu) + 1j*gamma) / (sigma *np.sqrt(2.0))
    return A * np.real(wofz(z))

def lorentzian(x, *p):
    A, p0, w, offset = p
    u = (p0 - x) / w
    return A * 1.0/(1.0+u**2) + offset


class LineFeature:
    """Each spectral line will be referenced as this object class.
    The line is automatically fit upon creation."""

    # reference available functions
    functions = [lorentzian,gaussian,voigt]
    
    def __init__(self, name, center, width, waves, flux,resample=None):
        self.name = '['+name+']'
        self.center = center
        self.line_text = '%s: %i' % (self.name,self.center)
        
        if width % 2 != 0:
            raise Exception('width must be an even number of pixels')
        self.width = width

        # Find the actual center of the feature
        self.center_idx, self.center_wave = self.find_center(waves,flux)

        # Calculate continuum
        self.waves,self.continuum,self.flux = self.find_continuum(waves,
                                                                  flux,
                                                                  width)
        # Integrate flux
        self.int_flux = integrate(self.waves,self.flux)
        
        # Fit function and integrate
        for function in LineFeature.functions:
            prefix = function.__name__

            waves, p = self.fit_function(function,width,resample)
            setattr(self, prefix+'_waves', waves)
            setattr(self, prefix+'_p', p)
            intflux = integrate(waves,function(waves,*p))
            setattr(self,prefix+'_flux',intflux)
            

        # Print results of fitting
        print '\t%s%i: Found feature at %.1f angstroms, cont = %.2e' % \
            (self.name,self.center,self.center_wave,self.continuum)
        print '\t\tFlux_int   = %.2e' % self.int_flux
        print '\t\tFlux_gauss = %.2e' % self.gaussian_flux
        print '\t\tFlux_loren = %.2e' % self.lorentzian_flux
        print '\t\tFlux_voigt = %.2e' % self.voigt_flux


        
    def find_center(self,waves,flux):
        """Find the center of the line by calculating local maxima"""
        local_maxima = argrelmax(flux)
        return find_nearest_element(waves[local_maxima],self.center,index=True)

    def find_continuum(self,waves,flux,width):
        """Calculate continuum around line center"""
        # Only want to fit data around the line
        idx = np.where(  (waves <= self.center_wave+width/2)
                        & (waves >= self.center_wave-width/2))[0]
        # Continuum should be outermost flux points
        cont = np.mean([flux[idx[0]],flux[idx[-1]]])
        return (waves[idx], cont,flux[idx]-cont)
        

    def fit_function(self,function,width,resample=None):
        """Apply specified function to data"""

        # Guess initial features [A,mu,sigma,gamma]
        p0 = [np.max(self.flux),self.center_wave,width/2.0,1.0]
        # Find best fitting params
        p, resid = curve_fit(function,self.waves,self.flux,p0=p0)

        ## If no resampling, just return the waves of the feature and params
        if resample == None:
            return (self.waves, p)

        ## Else, resample and return new waves
        else:
            line_x = np.linspace(self.center_wave-width/2.0,
                            self.center_wave+width/2.0,
                            num=resample)
            return (line_x, p)
            
    def plot(self, data=True, continuum=False,annotate=True):
        """Plot (but don't show) the feature"""
        # Reset color cycle for each feature
        plt.gca().set_color_cycle(None)

        # Plot original data
        if data:
            if continuum:
                plt.plot(self.waves,self.flux+self.continuum,'ko',markersize=8)
            else:
                plt.plot(self.waves,self.flux,'ko',markersize=8)
            
        # For each function, grab waves and params
        for function in LineFeature.functions:
            prefix = function.__name__
            x = getattr(self,prefix+'_waves')
            p = getattr(self,prefix+'_p')
            y = function(x,*p)
            # Plot the fitting function
            if continuum:
                plt.plot(x,y+self.continuum,linewidth=2)
                
            else:
                plt.plot(x,y,linewidth=2)
                            
        plt.legend(['data']+[x.__name__ for x in LineFeature.functions])
        plt.xlabel(r'$\lambda\,\,[\AA]$',fontsize=14)
        plt.ylabel(r'$F_{\lambda}\,\,[ergs/s/cm^2/\AA]$',fontsize=14)

        # Annotate line
        if annotate:
            if continuum:
                plt.text(self.center_wave,np.max(self.flux)+self.continuum,self.line_text,withdash=True,rotation='vertical',verticalalignment='bottom')
            else:
                plt.text(self.center_wave,np.max(self.flux),self.line_text,withdash=True,rotation='vertical',verticalalignment='bottom')
        #plt.show()


def integrate(x,y):
    """Integrate across y, excluding last edge"""
    dx = np.diff(x)
    # rip off last point of y
    y = y[:-1]
    return np.dot(y,dx)

        
def find_nearest_element(array,value,index=False):
    """Return closest value in array, or its index"""
    idx = np.abs(array-value).argmin()
    return (idx,array.flat[idx]) if index else array.flat[idx]


def get_1Dspec(filename):
    """Parse 1D spectrum from FITS file"""
    # Open FITS file and pull header
    ydata,header = pyfits.getdata(filename, header=True)

    # Construct wavelengths from header data
    waves = np.arange(0,len(ydata))*header['CDELT1'] + header['CRVAL1']
    waves -= 8.0  # to account for the broken offset

    # Create recarray, and return
    data = np.array(zip(waves,ydata),dtype=[('waves',np.float),
                                            ('flux',np.float)])
    return data



def main():
    """Client to LineFeature"""
    parser = argparse.ArgumentParser(description='Calculate line ratios of input file.')
    parser.add_argument('file',type=str,help='Input file for analysis')

    ## Parse command-line arguments
    args = parser.parse_args()

    ## Determine file type
    if os.path.splitext(args.file)[1] == '.txt':
        # If text file, simply read the columns
        data = np.genfromtxt(args.file, names=['waves', 'flux'])

    else:
        # Otherwise, read in FITS
        data = get_1Dspec(args.file)
        
    ## [OIII] lines
    print '[OIII]'    
    o5007 = LineFeature('OIII',5007,18,data['waves'],data['flux'],resample=100)
    o4959 = LineFeature('OIII',4959,16,data['waves'],data['flux'],resample=100)
    o4363 = LineFeature('OIII',4363,12,data['waves'],data['flux'],resample=100)
    oList = [o5007,o4959,o4363]

    oRatio = (o5007.gaussian_flux+o4959.gaussian_flux) / o4363.gaussian_flux
    print
    print '\t[j(5007)+j(4959)]/[j(4363)]: %.1f' % oRatio
    print

    ## [SII] lines
    print '[SII]'
    s6716 = LineFeature('SII',6716,12,data['waves'],data['flux'],resample=100)
    s6731 = LineFeature('SII',6731,12,data['waves'],data['flux'],resample=100)
    sList = [s6716,s6731]
    
    sRatio = s6716.gaussian_flux / s6731.gaussian_flux
    print
    print '\t[j(6716])/[j(6731)]: %.1f' % sRatio
    print

    ## Hydrogen
    print '[H]'
    hb = LineFeature('HBeta',4861,12,data['waves'],data['flux'],resample=100)
    ha = LineFeature('HAlpha',6563,16,data['waves'],data['flux'],resample=100)
    hList = [ha,hb]
    
    hRatio = ha.gaussian_flux / hb.gaussian_flux
    print
    print '\t[j(hA)]/[j(hB)]: %.1f' % hRatio
    print

    # Plot total spectrum
    plt.plot(data['waves'],data['flux'],'k',linewidth=1)
    
    # Plot all features and fits
    for feature in oList+sList+hList:
        feature.plot(data=False,continuum=True,annotate=False)

    # Show all lines
    plt.show()
    return 0

    
if __name__ == '__main__':
    main()
