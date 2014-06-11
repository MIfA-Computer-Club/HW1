"""


"""
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import os


def normexponent(val):
    """Return the exponent n such that 0 < val/10**n < 10."""
    n = np.log10(val)
    if n < 0:
        n = int(n) - 1
    else:
        n = int(n)
    return n


def load_data(filename):
    hdulist = fits.open(filename)
    hdu = hdulist[0]
    data, hdr = hdu.data, hdu.header
    return data, hdr


def get_wavelengths(data, hdr):
    wave0 = hdr['CRVAL1']
    dwave = hdr['CDELT1']
    wave = dwave*np.arange(data.size) + wave0
    return wave


def plot_spec(x, y, filename, logy=False, ymaxn=None):
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_axes([0.12, 0.12, 0.83, 0.83])

    if logy:
        y = np.log10(y)
        ymin, ymax = y.min(), y.max()
        if ymax < ymin:
            ymin, ymax = ymax, ymin
        dy = (ymax - ymin) * 0.05
        ymin, ymax = ymin - dy, ymax + dy
        ylabel = r'$\log_{10}[\mathrm{Flux / (\,erg \,s^{-1} \,cm^{-2} \,\AA^{-1})}]$'
    else:
        n = normexponent(y.max())
        y = y / 10**n
        ymin, ymax = 0, y.max()*1.05
        ylabel = (r'Flux ($10^{{{:d}}} '
                  '\mathrm{{\,erg \,s^{{-1}} \,cm^{{-2}} \,\AA^{{-1}}}}$)'
                  .format(n))

    ax.plot(x, y, 'k-')
    ax.set_ylim([ymin, ymax])

    if ymaxn:
        ax.yaxis.set_major_locator(plt.MaxNLocator(ymaxn))

    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel(ylabel)
    fig.savefig(filename)


def line_flux_simple(x, y, y0):
    """

    """
    # Rough continuum subtraction
    y = y - y0

    return None


def main():
    # 1)
    specfile = '~/Research/code/computer_club/HW1/hIIspec.fits'
    data, hdr = load_data(os.path.expanduser(specfile))

    # 2)
    wave = get_wavelengths(data, hdr)
    plotfile = '~/Research/code/computer_club/HW1/solutions/jake/spec_linear.pdf'
    plot_spec(wave, data, os.path.expanduser(plotfile))
    plotfile = '~/Research/code/computer_club/HW1/solutions/jake/spec_log.pdf'
    plot_spec(wave, data, os.path.expanduser(plotfile), logy=True, ymaxn=10)

    # 3)
    # Local continuum leves near 4959 and 5007

if __name__ == '__main__':
    main()
