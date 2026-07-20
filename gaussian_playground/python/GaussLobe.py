import UtilFunc as util
import BoysFunc as boys
import numpy as np
import scipy as sp

PI  = np.pi
N3D = 3

def gauss_prod(a,AV,b,BV):
    p  = a + b
    q  = (a*b) / (a+b)
    PV = (AV * (a/p)) + (BV * (b/p))
    sq_ab = util.calc_dist_sq(AV,BV)
    Kab   = np.exp(-q*sq_ab)
    return Kab,p,PV

def ovlap(a,AV,b,BV):
    p = (a+b)
    q = (a*b)/(a+b)
    sq_ab = util.calc_dist_sq(AV,BV)
    Kab   = np.exp(-q*sq_ab)
    Sab   = np.sqrt((PI/p)**(N3D)) * Kab
    return Sab

def ovlap_and_ke(a,AV,b,BV):
    p = (a+b)
    q = (a*b)/(a+b)
    sq_ab = util.calc_dist_sq(AV,BV)
    Kab   = np.exp(-q*sq_ab)
    Sab   = np.sqrt((PI/p)**(N3D)) * Kab
    z1    = 2.0e0*q*sq_ab
    z2    = (3.0e0-z1)*q
    Tab   = z2 * Sab
    return Sab,Tab

def vext(a,AV,b,BV,CV):
    p  = (a+b)
    q  = (a*b)/(a+b)
    PV = (AV * (a/p)) + (BV * (b/p))
    rsq_ab = util.calc_dist_sq(AV,BV)
    rsq_cp = util.calc_dist_sq(CV,PV)
    t1     = (2.0e0*PI/p)
    t2     = np.exp(-q*rsq_ab)
    t3     = p*rsq_cp
    t3     = boys.boys_func(t3)
    ans    = t1 * t2 * t3
    return ans

def vee(a,AV,b,BV,c,CV,d,DV):
    kab,pab,PVab = gauss_prod(a,AV,b,BV)
    kcd,pcd,PVcd = gauss_prod(c,CV,d,DV)
    rsq = util.calc_dist_sq(PVab,PVcd)
    z1 = pab*pcd*rsq/(pab+pcd)
    z2 = (pab*pcd)*np.sqrt(pab+pcd)
    z3 = 2.0e0*np.sqrt(PI**5)*boys.boys_func(z1)
    ans = z3*kab*kcd/z2
    return ans

def gem_ovlap(a1,AV1,b1,BV1,c2,CV2,d2,DV2,gem):
    kab,pab,PVab = gauss_prod(a1,AV1,b1,BV1)
    kcd,pcd,PVcd = gauss_prod(c2,CV2,d2,DV2)
    rsq = util.calc_dist_sq(PVab,PVcd)
    z1  = (pab*pcd)+(pab*gem)+(pcd*gem)
    z2  = np.sqrt((PI*PI/z1)**N3D)
    q   = pab*pcd*gem/z1
    k12 = z2*np.exp(-q*rsq)
    ans = k12*kab*kcd
    return ans
