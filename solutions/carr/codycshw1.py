##############################################################################################
                                   #Computer Club HW: 1
import math as m
import numpy as np
import scipy as s
from scipy.integrate import quad
import matplotlib.pyplot as plt 

##############################################################################################
                                      #open data file

x_list = []
y_list = []

with open('hii.ascii.txt', 'r') as f:
	for line in f:
		line = line.strip()
		line = line.split()
		x = float(line[0])
		x_list.append(x)
		y = float(line[1])
		y_list.append(y)

f.close()

##############################################################################################
                                        #plot

plt.plot(x_list,y_list, 'r')
plt.show()
