import scipy.special
import sys

#description:
#associated Legendre polynomials of first kind calculated via python

def legendrepmn(m,n,x):
    P_lmp, dP_lmp = scipy.special.lpmn(m,n,x)
    result = float(P_lmp[m][n])
    return result

m = int(sys.argv[1])
n = int(sys.argv[2])
x = float(sys.argv[3])
print legendrepmn(m,n,x)
