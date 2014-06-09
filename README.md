Astro Computer Club 2014
========================

Assignment 1, Due ~Friday, June 20
----------------------------------


I have attached a FITS file that is a 1D spectrum of an HII region that you will be analyzing.  This assignment is entirely ripped from one given by Terry in his ISM class (but restructured for your python education).  Your assignment has a few parts (and quickly read the whole email before starting).  This assignment should work on the physics network install of python, even if you have not yet completed the virtualenv setup in my previous email.

**This is a lot**.  Parse through it at your own speed.  Understanding the science of HII regions (and how/why you can use OIII/SII lines as tracers of temp/density) is much less important at this junction than just getting your coding skills developed.  (I'm sure every faculty member would completely disagree with the previous statement, but whatever).

Feel free to bother me about any of the steps below.  Most of the current grad students have done this before in Terry's class (but probably not using Python, so hopefully this will be a useful exercise for all), so please ask everyone questions.  This is supposed to be an exercise in learning how to code, learning how to ask constructive questions (of humans and Google), and learning how to learn.

------

0.5)  If you are new to python, consider spending some time with codeacademy to learn the basic syntax, or even just diddle with the language a bit before starting. Before you start the assignment, at the very least know how to run python in a terminal, import various packages, and make/run a script from a text file.

1)  Load the FITS file into python using the pyfits module (or astropy.io.fits).  Be sure to also load the header.  Google the pyfits documentation, and maybe look through the tutorial/quick start section.  (If you are struggling with the FITS shit, then just use the text file provided--look into parsing a whitespace-delimited text file by using the built-in ```open``` function, or perhaps better ```numpy.genfromtxt```)

2)  Plot the spectrum (check out the matplotlib.pyplot module---matplotlib is your best friend) using the appropriate wavelengths (you will have to get some information from the header to do this.  Check out [this](http://nbviewer.ipython.org/github/gabraganca/S4/blob/master/notebooks/load-spectrum-FITS.ipynb), lines 5-7).  Another helpful code snippet would be the following:
```python
waves = np.arange(0,len(ydata))*crdelt+crval
```

3)  Tell me the approximate electron temperature of the gas in the HII region using OIII line ratios (check out the attached images to see which lines I mean).
  * First, for the 4959 and 5007 OIII lines, integrate each emission feature to get the flux of the lines (but remember to subtract off the continuum flux [the flux just off the line]).  For the 4363 line (arrow pointing to it in [hiiOxygen.gif](../master/hiiOxygen.gif)), just choose the peak intensity of the line since it's so small/broad.  Terry messed up the wavelength calibration, so the whole spectrum is ~10 angstroms off (i.e. the 4959 emission line is at ~4967).  You should get flux values on the order of 10^-14 -> 10^-15 for the strong lines, and ~10^-16 for the 4363 line.
  * Use the formula on page 3 of this [PDF](http://www-astro.physics.ox.ac.uk/~pfr/C1_TT/C1_ISM_Lecture4.pdf) to calculate the electron temperature given the OIII line fluxes you calculated.

4)  Tell me the approximate electron density of the HII region using sulfur lines.
  * Again, find integrated fluxes of the emission features.  The two SII lines are at 6716 and 6731 angstroms.  Remember to subtract off the continuum flux, and remember that the wavelengths are shifted.
  * Refer to the graph from Osterbrock (available on page 11 of this [PDF](http://zuserver2.star.ucl.ac.uk/~msw/teaching/PHAS2521/notes_2.pdf)), to roughly determine the electron density based on your SII line ratio.  (The graph on page 10 of that same PDF should also reproduce the temp you calculated in 3b)

Bonus.  Rather than just integrate the flux in the emission lines, try to fit a Gaussian curve to the feature, resample to a higher resolution, and integrate that.  There are a ton of ways to do function fitting in python, but the easiest is perhaps ```scipy.optimize.curve_fit```, so look into that documentation.

Bonus+.  Try a Lorentzian profile fit, or even a Voigt profile, to better match the [spectral line shape](http://en.wikipedia.org/wiki/Spectral_line_shape).
