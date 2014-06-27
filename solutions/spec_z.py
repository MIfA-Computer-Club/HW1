  #! /usr/bin/env python
import numpy as np
import pyfits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

spectra=pyfits.open("HIIspec.fits")
data = spectra[0].data
hdr = spectra[0].header


pxscale = hdr['CDELT1']
start1 = hdr['CRVAL1']

wave=np.arange(0,len(data))*pxscale+ start1



def FindClosest(xlist,val):
    return np.argmin(np.abs(xlist - val))


wl=[4959., 5007.,4363.]
delta=15. #Angstroms
lstart=[e-delta+10. for e in wl]
lend=[e+delta+10. for e in wl]

start =[FindClosest(wave,s) for s in lstart]
stop = [FindClosest(wave,e) for e in lend]
ydat = [data[s:e] for s,e in zip(start,stop)]
xdat = [wave[s:e] for s,e in zip(start,stop)]


gaussfunc = lambda x,a,b,c,d: a * np.exp(-(x-b)**2/(2*c**2)) + d

p=[]
for i in range(len(ydat)):
    p.append([np.max(ydat[i]),xdat[i][len(xdat[i])//2],1.,0])
      
#def popt_pcov(xdat,ydat,p):
#	return popt_o3_I,pcov_o3_I = curve_fit(gaussfunc,xdat,ydat,p0)

#def flux():
#	return flux_line_o3_I=popt_o3_I[0]*popt_o3_I[2]*(2.*np.pi)**0.5

   
popt_o3_I,pcov_o3_I = curve_fit(gaussfunc,xdat[0],ydat[0],p0=p[0])
popt_o3_II,pcov_o3_II = curve_fit(gaussfunc,xdat[1],ydat[1],p0=p[1])
popt_un,pcov_un = curve_fit(gaussfunc,xdat[2],ydat[2],p0=p[2])


flux_line_o3_I=popt_o3_I[0]*popt_o3_I[2]*(2.*np.pi)**0.5
flux_line_o3_II=popt_o3_II[0]*popt_o3_II[2]*(2.*np.pi)**0.5
flux_line_un=popt_un[0]*popt_un[2]*(2.*np.pi)**0.5


print "OIII first : ",popt_o3_I, "FLUX---> ", flux_line_o3_I
print "OIII second : ",popt_o3_II, "FLUX---> ", flux_line_o3_II
print " 4363 : ",popt_un, "FLUX---> ", flux_line_un

plt.plot(wave,data,'b')
plt.xlim([4200,7500])
plt.plot(xdat[0],gaussfunc(xdat[0],*popt_o3_I),'r')
plt.plot(xdat[1],gaussfunc(xdat[1],*popt_o3_II),'r')
plt.plot(xdat[2],gaussfunc(xdat[2],*popt_o3_II),'r')
plt.show()


electronTemp=33000 / (np.log(0.14) - np.log(flux_line_un/(flux_line_o3_I+flux_line_o3_II)))

print "electron temprature is: ",electronTemp

