import mpmath
import sys

#description:
#calculation of the polylogarithm order n of x through python

n = int(sys.argv[1])
x = float(sys.argv[2])
print mpmath.polylog(n,x)
