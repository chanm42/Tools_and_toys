import l800_UtilFunc as util
import numpy as np
import scipy as sp

PI  = np.pi
N3D = 3

def boys_func(x):
    ans = 1.0e0
    if x > 0.0e0:
        t1 = (PI / (4.0e0*x))**0.5
        t2 = sp.special.erf(x**0.5)
        ans = t1 * t2
    return ans

