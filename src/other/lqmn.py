import scipy.special
import sys

#description:
#associated Legendre polynomials of second kind calculated via python

def legendreqmn(m,n,x):
    P_lmp, dP_lmp = scipy.special.lqmn(m,n,x)
    result = float(P_lmp[m][n])
    return result

m = int(sys.argv[1])
n = int(sys.argv[2])
x = float(sys.argv[3])
print legendreqmn(m,n,x)
