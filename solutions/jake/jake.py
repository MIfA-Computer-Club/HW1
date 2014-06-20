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


def line_flux_simple(x, y, x1, x2, y0=None):
    """Integrate y as a function of x between x1 and x2

    Simple trapezoid integrator with an option for constant continuum
    subtraction.

    Parameters
    ----------
    x, y : array
        x and y data.
    x1, x2 : float
        Integration limits.
    y0 : float, optional
        If specified, subtract off this amount from the y values before
        integration. Default is None.

    Returns
    -------
    float

    """
    # Only want the data within the integration limits
    idx = (x1 <= x) & (x <= x2)
    x, y = x[idx], y[idx]

    # Rough continuum subtraction
    y = y if y0 is None else y - y0

    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    area = np.sum(dx*y[:-1] + 0.5*dx*dy)

    return area


def electron_temp(I4363, I4959, I5007):
    """Return the electron temperature (K) for the given line strengths."""
    return 33000 / (np.log(0.14) - np.log(I4363/(I4959+I5007)))


def electron_density(I6716, I6731):
    """Return the electron density (cm3) from the given line strengths.

    The ratio vs. density data were lifted from the Osterbrock plot using
    PlotDigitizer.

    """
    density_data = np.array([
        9.866987, 19.478735, 36.90968, 72.14769, 122.13107, 226.92259, 376.6114,
        612.408, 965.6635, 1538.101, 2687.6423, 4993.7065, 9865.964, 19285.09,
        37684.582, 75936.086, 100297.75
        ])
    ratio_data = np.array([
        1.4327427, 1.4224162, 1.3981028, 1.3626771, 1.3047774, 1.2219033,
        1.1138053, 0.99453276, 0.8668605, 0.75313777, 0.64510214,
        0.5622281, 0.5045155, 0.4690898, 0.453176, 0.44844922, 0.44321114,
        ])

    # reverse so that ratio_data is increasing
    ratio_data = ratio_data[::-1]
    density_data = density_data[::-1]

    return np.interp(I6716/I6731, ratio_data, density_data)


def main():
    # 1)
    specfile = '~/Research/code/computer_club/HW1/hIIspec.fits'
    data, hdr = load_data(os.path.expanduser(specfile))

    # 2)
    wave = get_wavelengths(data, hdr)
    plotfile = '~/Research/code/computer_club/HW1/solutions/jake/spec.pdf'
    plot_spec(wave, data, os.path.expanduser(plotfile))

    # 3)
    I4959 = line_flux_simple(wave, data, 4950, 4990, 2.8e-16)
    I5007 = line_flux_simple(wave, data, 4990, 5040, 2.8e-16)
    I4363 = line_flux_simple(wave, data, 4360, 4380, 3.0e-16)
    Te = electron_temp(I4363, I4959, I5007)
    print 'electron temp (K) = {:.2e}'.format(Te)  # 1.31e4

    # 4)
    I6716 = line_flux_simple(wave, data, 6715, 6734, 1.8e-16)
    I6731 = line_flux_simple(wave, data, 6734, 6750, 1.8e-16)
    ne = electron_density(I6716, I6731)
    print 'electron density (cm3) = {:.2e}'.format(ne)  # 3.27e1


if __name__ == '__main__':
    main()
