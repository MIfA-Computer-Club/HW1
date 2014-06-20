#Computer Club Assignment 1

#6/20/16

#Part 2: Plot spectrum

import pyfits
from pylab import *
import numpy as np
import matplotlib
import os
plt.rc('text', usetex=True)


spec=pyfits.open('hIIspec.fits')

flux=spec[0].data

ref_pixel=spec[0].header['CRPIX1']
coord_ref_pixel=spec[0].header['CRVAL1']
wave_pixel=spec[0].header['CDELT1']

def get_wstart(ref,wave_ref,wave_per_pixel):

    return wave_ref - ((ref-1)*wave_per_pixel)

wstart=get_wstart(ref_pixel,coord_ref_pixel,wave_pixel)

x=np.linspace(wstart,wstart+wave_pixel*len(flux),len(flux))   
plt.plot(x,flux)
xlabel(r'$\mathrm{Wavelength}~\AA$', fontsize=20)
ylabel(r'$\mathrm{Flux}$', fontsize=20)
#xlim(4800,5030)
plt.show()

#Part 3: Get peak flux from 4363, 4959 and 5007 OIII lines, use to find electron temperature. 4363 is small enough to just use peak value.
terryerror=8
hw=10
#4363 peak:
continuum=average(flux)
r4363=[]
for i in range(int((4363+terryerror-hw-wstart)/1.2),int((4363+terryerror+hw-wstart)/1.2)):
    r4363.append(flux[i]-continuum)

flux4363=max(r4363)*1.2
#janky test method of integration - pseudo-midpoint method, multiply each flux value by width (1.2 A)
#4959 integrated
r4959=[]
for i in range(int((4959+terryerror-hw-wstart)/1.2),int((4959+terryerror+hw-wstart)/1.2)):
    r4959.append(flux[i])

flux4959=sum(r4959-continuum)*1.2

#5007 integrated
r5007=[]
for i in range(int((5007+terryerror-hw-wstart)/1.2),int((5007+terryerror+hw-wstart)/1.2)):
    r5007.append(flux[i])

flux5007=sum(r5007-continuum)*1.2

I=flux4363/(flux4959+flux5007)

Temp=-33000./log(I/.14)
