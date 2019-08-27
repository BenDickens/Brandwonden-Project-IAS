import math as m
import numpy as np

a = [(1,2), (3,4), (2,8)]
b = (3,3)
e = []
for x in a:
	c = b[0]-x[0]
	d = b[1]-x[1]
	e.append((c,d))

h = 10000
for y in e:
	p = y[0]
	q = y[1]
	if(m.hypot(p,q)<h):
		h = m.hypot(p,q)
	else

