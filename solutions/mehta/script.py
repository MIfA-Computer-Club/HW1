import pyfits
import numpy as np
import scipy.optimize
import scipy.integrate
import matplotlib.pyplot as plt
from scipy.special import wofz

class Line:
    """
    Class for spectral lines
    Fits the line when first called
    Provide a name, line center, wavelength and flux arrays
    Array size doesn't matter as long as it contains line center
    Class will have line name, fitted center, norm, FWHM and flux after init
    """
    def __init__(self,name,center,waves,spec):

        window, cdelt = 20., waves[1]-waves[0]
        spec_cut = spec[ ((center-0.5*window) <= waves) & (waves <= (center+0.5*window)) ]
        [center] = np.where(spec == max(spec_cut))
        spec_cut = spec[center-int(0.5*window/cdelt):center+int(0.5*window/cdelt)]
        waves_cut = waves[center-int(0.5*window/cdelt):center+int(0.5*window/cdelt)]
        popt, pcov = scipy.optimize.curve_fit(voigt_fwhm, waves_cut, spec_cut, p0=(np.mean(spec_cut),np.mean(waves_cut),1.0,1.0))
        self.name, [self.A,self.x0,self.fwhm,self.gamma], [self.dA,self.dx0,self.dfwhm,self.dgamma] = name, popt, np.sqrt(np.diag(pcov))
        self.flux, self.dflux = scipy.integrate.quad(voigt_fwhm, min(waves_cut)-50, max(waves_cut)+50, args=(self.A,self.x0,self.fwhm,self.gamma))

def gauss_fwhm(x, A, x0, fwhm):
    """
    Gaussian function using FHWM
    """
    sig = fwhm/2./np.sqrt(2*np.log(2)) #Convert FWHM to sigma
    return A*np.exp(-(x-x0)**2/2./sig**2)
    
def voigt_fwhm(x, A, mu, fwhm, gamma):
    """
    Voigt function using FWHM
    """
    sigma = fwhm/2./np.sqrt(2*np.log(2)) #Convert FWHM to sigma
    z = ((x-mu) + 1j*gamma) / (sigma *np.sqrt(2.0))
    return A * np.real(wofz(z))
    
def get_data(fname):
    """
    Get data from FITS spectrum
    """
    data = pyfits.open(fname)
    hdr = data[0].header
    spec = data[0].data
    crpix1, cdelt1, crval1 = hdr['CRPIX1'], hdr['CDELT1'], hdr['CRVAL1']
    data.close()

    waves = ((np.array(range(spec.size))+crpix1) * cdelt1) + crval1
    waves = waves - 10 ### Correction for Calibration Error
    
    return waves, spec
    
def continuum_fit(waves,spec):
    """
    Continuum Fitting
    """
    wave_cont, spec_cont = waves, spec
    wave_cont, spec_cont = wave_cont[spec_cont < np.std(spec_cont)], spec_cont[spec_cont < np.std(spec_cont)] # Remove parts that deviate more than 1 sigma

    #Smooth spectrum by box-car
    w = np.array([0,0,1,1,1,0,0]) # 3 px Box
    spec_cont = np.convolve(w,spec_cont/w.sum(),'same') # Convolve spectrum
    cont_fit = np.polyfit(wave_cont,spec_cont,4) # Fit continuum
    cont_func = np.poly1d(cont_fit) # Create continuum function
    return cont_func
    
def plot():
    """
    Main plotting function
    """
    waves, spec = get_data("../../hIIspec.fits") # Get data

    fig, ax = plt.subplots(1,1,figsize=(20,8),dpi=100) # Initialize figure
    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel(r'Flux [$erg \cdot\ s^{-1} \cdot\ cm^{-2}$]')
    ax.set_title('Spectrum')
        
    ax.plot(waves,spec,c='k',lw=1.0,alpha=0.5) # Plot the spectrum
    continuum = continuum_fit(waves,spec) # Fit continuum
    spec = spec - continuum(waves) # Subtract continuum off

    # Fit the 5 lines in question
    # Line variables are made global so they may be used
    # by other funcs without explicitly passing as args
    global o4363,o4959,o5007,s6716,s6731
    o4363 = Line('[OIII]4363', 4363., waves, spec)
    o4959 = Line('[OIII]4959', 4959., waves, spec)
    o5007 = Line('[OIII]5007', 5007., waves, spec)
    s6716 = Line('[SII]6716', 6716., waves, spec)
    s6731 = Line('[SII]6731', 6731., waves, spec)
    lines = [o4363,o4959,o5007,s6716,s6731]
    
    # Setup a colormap for the lines
    colormap = plt.cm.gist_rainbow
    cindex = np.linspace(0.1,0.9,len(lines))

    # Overplot the fitted lines on the spec
    for i,line in enumerate(lines):
    
        plot_waves = np.linspace(line.x0-25.,line.x0+25.,100) # High resolution wavelength range
        ax.plot(plot_waves, voigt_fwhm(plot_waves,line.A,line.x0,line.fwhm,line.gamma)+continuum(plot_waves), c=colormap(cindex[i]), label=line.name)
        
        # Output the parameters for each line
        print("Line:\t%s" % (line.name))
        print(u"Fit:\tCenter\t= %.2f \u00B1 %.3f [Angs]" % (line.x0,line.dx0))
        print(u"\tSigma\t= %.4f \u00B1 %.4f [Angs]" % (line.fwhm,line.dfwhm))
        print(u"\tGamma\t= %.4f \u00B1 %.4f [Angs]" % (line.gamma,line.dgamma))
        print(u"Flux:\t%.2e [ergs/s/cm^-2]\n" % (line.flux))

    ax.legend() # Enable the legend
    fig.savefig('spectrum.png') # Save plot

def calc_T(o4363,o4959,o5007): 
    """
    Calculate the Temp. using OIII line ratios
    """
    o3_ratio = o4363/(o4959+o5007)
    T = -33000./np.log(o3_ratio/0.14)
    print("[OIII]4363/[OIII]4959+5007 = %.4f" % o3_ratio)
    print("Temperature: %.2f K \n" % T)

def calc_ne(s6716,s6731):
    """
    Calculate the electron density using SII line ratios
    No mathematical expr. for now
    Using direct value from graph
    """
    s2_ratio = s6716/s6731
    n_e = 40. # from graph
    print("[SII]6716/[SII]6731 = %.4f" % s2_ratio)
    print("e Density: %.2f cm^-3\n" % n_e)


plot()
calc_T(o4363.flux,o4959.flux,o5007.flux)
calc_ne(s6716.flux,s6731.flux)
