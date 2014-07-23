import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

def main():

    xdata = np.linspace(0, 4, 50)
    y = func(xdata, 2.5, 1.3, 0.5)
    ydata = y + 0.2 * np.random.normal(size=len(xdata))

    popt, pcov = curve_fit(func, xdata, ydata)
    
    plt.plot(xdata,ydata)
    plt.plot(xdata,func(xdata, popt[0], popt[1], popt[2]))
    plt.show()
    
main()
