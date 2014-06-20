##############################################################################################
                                   #Computer Club HW: 1
import math as m
"""JS: I wouldn't bother importing the math module if you're just going to
import numpy anyway, which does all that math can do and more.

"""
import numpy as np
import scipy as s
from scipy.integrate import quad
import matplotlib.pyplot as plt 

##############################################################################################
                                      #open data file

x_list = []
y_list = []
"""JS: For what it's worth, there's a strong preference in the python
community for spaces over tabs, usually 4 spaces per indentation level.

"""
with open('hii.ascii.txt', 'r') as f:
	for line in f:
		line = line.strip()
                """JS: There's little reason to strip off whitespace in
                this case. The float conversion will take care of any
                whitespace remaining after splitting the line.

                """
		line = line.split()
		x = float(line[0])
		x_list.append(x)
		y = float(line[1])
		y_list.append(y)

f.close()
"""JS: You don't have to explicity close files when you use the `with`
statement.

"""

##############################################################################################
                                        #plot

plt.plot(x_list,y_list, 'r')
plt.show()
