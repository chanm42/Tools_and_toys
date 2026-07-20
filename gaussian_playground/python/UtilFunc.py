import numpy as np
PI  = np.pi
N3D = 3

def calc_dist_sq(xvec,yvec):
    ndim   = len(xvec)
    distsq = 0.0
    for i in range(ndim):
        distsq = distsq + ((xvec[i] - yvec[i])**2)
    return distsq

def calc_dist(xvec,yvec):
    distsq = calc_dist_sq(xvec,yvec)
    dist   = distsq**0.5
    return dist