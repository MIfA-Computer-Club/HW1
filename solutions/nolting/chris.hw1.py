import numpy
import pyfits
import matplotlib.pyplot as plt
import scipy.integrate as i
from scipy.optimize import curve_fit

info = 0



def gauss(x,a,b,c,d):
    return a*numpy.exp(-(x-b)**2/(2*c**2))+d



def line(wave,dat,minWave, maxWave,lineNum):
	
	array1 = numpy.where((wave>minWave)&(wave<maxWave))
	wave1 = wave[array1]
	dat1 = dat[array1]
	cont1 = 0.5*(dat1[0]+dat1[len(dat1)-1])
	line1 = i.simps(dat1-cont1,wave1)
	
	plt.plot(wave1, dat1-cont1)
	plt.xlabel('Wavelength')
	plt.ylabel('Intensity')
	plt.savefig(lineNum)
	plt.clf()
	
	p0 = [numpy.amax(dat1),(maxWave+minWave)/2.0,(maxWave-minWave)/2.35,cont1]
	wave2 = numpy.linspace(minWave,maxWave,100)
	
	try: 
		param , cov = curve_fit(gauss,wave1,dat1,p0,maxfev=10000)
	except RuntimeError:
		param = [1,1,1,1]
		print("Error - curve_fit failed"), lineNum
    
	print 'FINAL PARAM', param
	[a1, b1, c1, d1] = param	
	gaussfit = gauss(wave2,a1,b1,c1,d1)
	
	plt.plot(wave1,dat1)
	plt.plot(wave2,gaussfit)
	plt.xlabel('Wavelength')
	plt.ylabel('Intensity')
	gaussname = lineNum + '_gaussfit'
	plt.savefig(gaussname)
	plt.clf()
	
	line2 = i.simps(gaussfit-cont1,wave2)
	
	return line1, line2

def main():

	h2spec = pyfits.open('hIIspec.fits')
	head = h2spec[0].header
	dat = h2spec[0].data

	if info == 1:
		h2spec.info()
		print repr(head)

	refPix = h2spec[0].header['CRPIX1']
	refPixCoord = h2spec[0].header['CRVAL1']
	delta = h2spec[0].header['CDELT1']

	wave = numpy.arange(refPixCoord,refPixCoord + delta* len(dat),delta)

#	##use this to find the range over which to integrate the relavant lines
#	plt.plot(wave,dat) 
#	plt.show()
#	
	Ox4363, Ox4363g = line(wave,dat,4360,4380,'Ox4363')
	Ox4959, Ox4959g = line(wave,dat,4960,4974,'Ox4959')
	Ox5007, Ox5007g = line(wave,dat,5008,5021,'Ox5007')
	
	OxRatio = Ox4363/(Ox5007+Ox4959)
	OxGRatio = Ox4363g/(Ox5007g+Ox4959g)
	print 'Oxygen ratio is ', OxRatio #~0.01078 gives me T_e ~ 10^5 K
	print 'Gaussian fit oxygen ratio is ', OxGRatio
	
	S6716, S6716g = line(wave,dat,6720,6733,'S6716')
	S6731, S6731g = line(wave,dat,6735,6748,'S6731')
	SRatio = S6716/S6731
	SGRatio = S6716g/S6731g
	print 'Silicon ratio is ', SRatio #~1.383 gives me n_e ~ 50 cm^-3
	print 'Gaussian fit silicon ratio is ', SGRatio, 'DO NOT TRUST - BROKEN'
	
	h2spec.close()

main()
