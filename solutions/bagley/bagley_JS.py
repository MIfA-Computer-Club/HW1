#! /usr/bin/env python
import argparse
import pyfits
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.special import wofz


def get_peak(wave, spec, wavemin, wavemax):
    """Create arrays for which wave[wavemin:wavemax] and 
    spec[wavemin:wavemax]. Find the index of the line peak."""
    xx = wave[np.where( (wave > wavemin) & (wave < wavemax) )]
    yy = spec[np.where( (wave > wavemin) & (wave < wavemax) )]
    """JS: more compact, if you prefer::

      xx = wave[(wave > wavemin) & (wave < wavemax)]

    The stuff in the brackets returns a boolean array, which indexes
    `wave`.

    """
    peak = np.argmax(yy)
    return {'peak':peak, 'xx':xx, 'yy':yy}
 

def get_continuum(yy, peak_idx, width):
    """Calculate the continuum as the average of the flux 
    at +/- width/2 from line center."""
    return np.mean( [yy[peak_idx - width/2.], yy[peak_idx + width/2.]] )

   
def integrate_flux(xx, yy):
    """Integrate line using Simpson's rule"""
    # If there are an odd number of intervals, take average of:
    #   1) first N-2 intervals with trapezoidal rule on last interval
    #   2) last N-2 intervals with trapezoidal rule on first interval
    return integrate.simps(yy, xx, even='avg')    


def fit_gauss(xx, a, mu, sigma, d):
    """Fit a Gaussian to the line."""
    return a * np.exp( -(xx-mu)**2 / (2.*sigma**2)) + d


def fit_lorentz(xx, a, p0, w, d):
    """ Fit a Lorentzian to the line."""
    r = (p0 - xx) / (w/2.)
    return a * 1./(1. + r**2) + d


def fit_voigt(xx, a, nu, sigma, gamma):
    """Fit a Voigt function to the line."""
    # The Voigt function:
    #   V(x;sigma,gamma) = Re[w(z)] / (sigma * sqrt(2 pi))
    #   where w(z) = exp(-z**2) * erfc(-i * z), the Faddeeva function
    #   evaluated for z = (x + i*gamma) / (sigma * sqrt(2))
    #   Use wofz() from scipy
    z = ((xx - nu) + 1j*gamma) / (sigma * np.sqrt(2.))
    return a * np.real(wofz(z)) / (sigma * np.sqrt(2. * np.pi))


def use_func(function, xx, yy, continuum, peak_idx, fwhm, res):
    """Fit the line with a given function (Gaussian, Lorentzia, Voigt"""
    # initial guess for parameters
    p0 = [np.max(yy), xx[peak_idx], fwhm/2., 1.]
    
    # fit function
    popt,pcov = curve_fit(function, xx, yy-continuum, p0=p0)

    # resample 
    step_size = (xx[-1:] - xx[0]) / res
    newx = np.arange(xx[0], xx[-1:], step_size)

    # integrate at new resolution
    flux = integrate_flux(newx, function(newx, popt[0], popt[1],
                                               popt[2], popt[3]) )

    # plot lines to check
    plt.plot(xx,yy)
    plt.plot(newx, function(newx, *popt)+continuum, 'r')
    plt.plot([xx[0], xx[-1:]], [continuum, continuum], 'k')
    plt.show()
    return flux

def fit_line(wave, spec, wavemin, wavemax, scale, res):
    """Get flux of line by simple integration, and fitting
    Gaussian, Lorentzian and Voigt profiles."""
    # Begin by fitting Gaussian to get line width (for finding continuum)
    # initial guess for parameters
    result = get_peak(wave, spec, wavemin, wavemax)
    peak_idx = result['peak']
    xx = result['xx']
    yy = result['yy']
    p0 = [np.max(yy), xx[peak_idx], 10., 1.]
    popt,pcov = curve_fit(fit_gauss, xx, yy, p0=p0)
    # FWHM ~ 2.355 * sigma
    fwhm = 2.355 * popt[2]
    # width is wide enough to reach continuum on either side of line
    width = scale * fwhm

    # now find continuum
    continuum = get_continuum(yy, peak_idx, width)

    # get line fluxes
    # basic integration
    int_flux = integrate_flux(xx, yy-continuum)
    
    # fit with gaussian
    gauss_flux = use_func(fit_gauss, xx, yy, continuum, peak_idx, fwhm, res)

    # fit with lorentzian
    lorentz_flux = use_func(fit_lorentz, xx, yy, continuum, peak_idx, fwhm, 
                                                                        res)

    # fit with voigt function
    voigt_flux = use_func(fit_voigt, xx, yy, continuum, peak_idx, fwhm, res)

    return {'integrated':int_flux, 'gaussian':gauss_flux, 
            'lorentzian':lorentz_flux, 'voigt':voigt_flux}
    
def line_ratios(wave, spec):
    # OIII
    f4959 = fit_line(wave, spec, 4930., 5000., 6, 100.)
    f5007 = fit_line(wave, spec, 4980., 5050., 6, 100.)
    f4363 = fit_line(wave, spec, 4350., 4390., 2, 100.)
    
    print '\nOIII'
    print '  [OIII]4959 fluxes: '
    print '    Simple integration = %.3e' % f4959['integrated']
    print '          Gaussian fit = %.3e' % f4959['gaussian']
    print '        Lorentzian fit = %.3e' % f4959['lorentzian']
    print '             Voigt fit = %.3e' % f4959['voigt']
    print '  [OIII]5007 fluxes: '
    print '    Simple integration = %.3e' % f5007['integrated']
    print '          Gaussian fit = %.3e' % f5007['gaussian']
    print '        Lorentzian fit = %.3e' % f5007['lorentzian']
    print '             Voigt fit = %.3e' % f5007['voigt']
    print '  [OIII]4363 fluxes: '
    print '    Simple integration = %.3e' % f4363['integrated']
    print '          Gaussian fit = %.3e' % f4363['gaussian']
    print '        Lorentzian fit = %.3e' % f4363['lorentzian']
    print '             Voigt fit = %.3e' % f4363['voigt']
    """JS: The above is certainly readable, but there are many different
    ways to do this. You can also do something like this::

      print (
          '\nOIII\n'
          '  [OIII]4959 fluxes:\n'
          '    Simple integration = %.3e\n'
          '          Gaussian fit = %.3e\n'
          '        Lorentzian fit = %.3e\n'
          '             Voigt fit = %.3e\n'
          '  [OIII]5007 fluxes:\n'
          '    Simple integration = %.3e\n'
          '          Gaussian fit = %.3e\n'
          '        Lorentzian fit = %.3e\n'
          '             Voigt fit = %.3e\n'
          '  [OIII]4363 fluxes:\n'
          '    Simple integration = %.3e\n'
          '          Gaussian fit = %.3e\n'
          '        Lorentzian fit = %.3e\n'
          '             Voigt fit = %.3e'
          ) % (f4959['integrated'], f4959['gaussian'],
               f4959['lorentzian'], f4959['voigt'],
               f5007['integrated'], f5007['gaussian'],
               f5007['lorentzian'], f5007['voigt'],
               f4363['integrated'], f4363['gaussian'],
               f4363['lorentzian'], f4363['voigt'])

    or this (using "new-style" formatting)::

      message = (
          '\nOIII',
          '  [OIII]4959 fluxes: ',
          '    Simple integration = {0:.3e}'.format(f4959['integrated']),
          '          Gaussian fit = {0:.3e}'.format(f4959['gaussian']),
          '        Lorentzian fit = {0:.3e}'.format(f4959['lorentzian']),
          '             Voigt fit = {0:.3e}'.format(f4959['voigt']),
          '  [OIII]5007 fluxes: ',
          '    Simple integration = {0:.3e}'.format(f5007['integrated']),
          '          Gaussian fit = {0:.3e}'.format(f5007['gaussian']),
          '        Lorentzian fit = {0:.3e}'.format(f5007['lorentzian']),
          '             Voigt fit = {0:.3e}'.format(f5007['voigt']),
          '  [OIII]4363 fluxes: ',
          '    Simple integration = {0:.3e}'.format(f4363['integrated']),
          '          Gaussian fit = {0:.3e}'.format(f4363['gaussian']),
          '        Lorentzian fit = {0:.3e}'.format(f4363['lorentzian']),
          '             Voigt fit = {0:.3e}'.format(f4363['voigt'])
          )  # could also do all the formatting here at the end as in the last example
      print '\n'.join(message)

    or my personal favorite::

      message = ('\nOIII\n'
                 '  [OIII]4959 fluxes:\n{0:s}'
                 '  [OIII]5007 fluxes:\n{1:s}'
                 '  [OIII]4363 fluxes:\n{2:s}')
      submsg = ('    Simple integration = {0:.3e}\n'
                '          Gaussian fit = {1:.3e}\n'
                '        Lorentzian fit = {2:.3e}\n'
                '             Voigt fit = {3:.3e}\n')
      key_list = ['integrated', 'gaussian', 'lorentzian', 'voigt']
      flist = [f4959, f5007, f4363]
      submsg_list = [submsg.format(*[f[key] for key in key_list]) for f in flist]
      print  message.format(*submsg_list)

    etc.

    """
    # use voigt integration for ratios
    oiii_ratio = (f4959['voigt'] + f5007['voigt']) / f4363['voigt']
    print '\n  j_4363 / (j_4959 + j_5007) = %.2f' % oiii_ratio
    # get electron temperature
    temp = -33000. / np.log((1./oiii_ratio) / 0.14)
    print '  Electron temperature = %.0f' % temp
    
    # SII
    f6716 = fit_line(wave, spec, 6700., 6740., 2.5, 100.)
    f6731 = fit_line(wave, spec, 6730., 6770., 2, 100.)

    print '\n\nSII'
    print '  [SII]6716 fluxes: '
    print '    Simple integration = %.3e' % f6716['integrated']
    print '          Gaussian fit = %.3e' % f6716['gaussian']
    print '        Lorentzian fit = %.3e' % f6716['lorentzian']
    print '             Voigt fit = %.3e' % f6716['voigt']
    print '  [SII]6731 fluxes: '
    print '    Simple integration = %.3e' % f6731['integrated']
    print '          Gaussian fit = %.3e' % f6731['gaussian']
    print '        Lorentzian fit = %.3e' % f6731['lorentzian']
    print '             Voigt fit = %.3e' % f6731['voigt']
    sii_ratio = f6716['gaussian'] / f6731['gaussian']
    print '\n  j_6716 / j_6731 = %.2f' % sii_ratio
    


def main():
    parser = argparse.ArgumentParser(
             description='Find line ratios of input spectrum')
    parser.add_argument('filename', type=str, nargs=1,
        help='Input spectrum file, either .fits or .txt')
    args = parser.parse_args()

    """JS: the `os` module is great for things like this::

      if os.path.splitext(filename[0])[1] == '.txt':
          # etc.

    This way you don't have to worry about the character length of the
    extension.

    """
    # is file .fits or .txt?
    if args.filename[0][-4:] == '.txt':
        # read columns
        data = np.genfromtxt(args.filename[0], names=['wave', 'spec'], 
                                            dtype=(float,float))
        """JS: python style police - wrap long function calls like this::

          data = np.genfromtxt(args.filename[0], names=['wave', 'spec'],
                               dtype=(float,float))

        or this::

          data = np.genfromtxt(
              args.filename[0], names=['wave', 'spec'], dtype=(float,float))
          #   ^ one or two indents to separate from what comes next

        """
        wave = data['wave']
        spec = data['spec']

    elif args.filename[0][-4:] == 'fits':
        # read in fits file and header
        spec,hdr = pyfits.getdata(args.filename[0], header=True)

        # get wavelength 
        wave = np.arange(0, len(spec)) * hdr['cdelt1'] + hdr['crval1']
    else:
        print 'File must be either .fits or .txt'
        exit()

    # plot the spctrum
    plt.plot(wave, spec, 'b')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel(r'flux (ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$ sr$^{-1}$) ')
    plt.ticklabel_format(useOffset=False)
    plt.savefig('bagley_spectrum.png')
    """JS: `savefig` has a dpi keyword for rasterized elements. Personally,
    I find the default (80, I think) makes for grainy text. Bumping it up to
    100 or more helps a lot, e.g., ::

      plt.savefig('bagley_spectrum.png', dpi=150)

    Of course, you could also just use a vector format (svg, pdf) instead.
    Probably ideal for a simple line plot like this.

    """
    plt.show()
    plt.close()

    # get line ratios
    line_ratios(wave, spec)


if __name__ == '__main__':
    main()
