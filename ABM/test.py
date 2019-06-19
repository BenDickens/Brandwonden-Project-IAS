import numpy as np
from matplotlib import pyplot
radius = 5
x0 = 5; a = radius  # x center, half width
y0 = 5; b = radius  # y center, half height
x = np.linspace(0, 10, 100)  # x values of interest
y = np.linspace(0, 10, 10)[:,None]  # y values of interest, as a "column" array
ellipse = ((x-x0)/a)**2 + ((y-y0)/b)**2 <= 1
print(ellipse)