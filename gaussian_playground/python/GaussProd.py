import UtilFunc as util
import BoysFunc as boys
import GaussLobe as gauss
import numpy as np
import scipy as sp

PI = np.pi
N3D = 3

def gauss_prod(a,AV,b,BV):
    p  = a + b
    q  = (a*b) / (a+b)
    PV = (AV * (a/p)) + (BV * (b/p))
    sq_ab = util.calc_dist_sq(AV,BV)
    Kab   = np.exp(-q*sq_ab)
    return Kab,p,PV

def gauss_triple_prod_1d(a,A,b,B,c,C):
    gam   = a+b  
    P     = (a*A + b*B) / gam
    q1    = (a*b) / gam
    K1    = np.exp(-q1*(A-B)**2)
    delta = gam + c 
    Q     = (gam*P + c*C) / delta
    q2    = (gam*c) / delta
    K2    = np.exp(-q2*(P-C)**2)
    ans   = K1*K2*(np.sqrt(PI / (a+b+c)))
    return ans

def gauss_triple_prod_1d_xy(a,A,b,B,c):
    gam1 = a+b 
    q    = (a*b) / gam1
    ans  = PI * np.sqrt((1.0e0 / (a*c + b*c)) * np.exp(-q * (A-B)**2))
    return ans

def gauss_triple_prod_3d_xy(a,AV,b,BV,c):
    ndim = len(AV)
    assert ndim == len(BV), ">>Error: Dimension mismatch in gauss_triple_prod_3d_xy"
    ans = 1.0e0
    for idim in range(ndim):
        centA = AV[idim]
        centB = BV[idim]
        ans   = ans * gauss_triple_prod_1d_xy(a,centA,b,centB,c)
    return ans

def gauss_triple_prod_1d_xy_zeros(a,b,c):
    gam = a+b+c
    KAB = np.sqrt(PI / gam)
    new_alp = (c + c**2) / gam
    return KAB,new_alp

def gauss_triple_prod_1d_xy_new(a,A,b,B,c):
    gam1 = a + b
    gam2 = a + b + c
    new_alp  = (gam1 * c) / gam2
    q        = a*b / gam1
    new_cent = (a*A + b*B) / gam1
    coef     = np.sqrt(PI / gam2)
    KAB      = coef * np.exp(-q * (A-B)**2)
    return new_alp,new_cent,KAB

def vee(a1,AV1,b1,BV1,c2,CV2,d2,DV2):
    kab,pab,PVab = gauss_prod(a1,AV1,b1,BV1)
    kcd,pcd,PVcd = gauss_prod(c2,CV2,d2,DV2)
    rsq = util.calc_dist_sq(PVab,PVcd)
    z1  = pab*pcd*rsq/(pab+pcd)
    z2  = (pab*pcd)*np.sqrt(pab+pcd)
    z3  = 2.0e0*np.sqrt(PI ** 5)*boys.boys_func(z1)
    ans = z3*kab*kcd/z2
    return ans

def full_contraction_zeros_1d(a1,b2,c3,d4,e13,i1,j2,k3,l4,m24):
    alp_13,K13 = gauss_triple_prod_1d_xy_zeros(c3,k3,e13)
    alp_24,K24 = gauss_triple_prod_1d_xy_zeros(d4,l4,m24)
    # now we have <a1,i1,alp_13,||b2,j2,alp_24>
    # contract the (a1,i1) gaussian and (b2,j2) gaussian using 0 for the center
    K11,alp_11,cent_11 = gauss_prod(a1,0,i1,0)
    K22,alp_22,cent_22 = gauss_prod(b2,0,j2,0)
    # now calculate the vee using <alp_11,alp_13|r12|alp_22,alp_24> and multiply by the coefficients K11*K13*K22*K24
    kern       = vee(alp_11,0,alp_13,0,alp_22,0,alp_24,0)
    ans        = K13 * K11 * K24 * K22 * kern
    return ans

def full_contraction_1d(a1,A,b2,B,c3,C,d4,D,e13,i1,I,j2,J,k3,K,l4,L,m24):
    alp_13,cent_13,K13 = gauss_triple_prod_1d_xy(c3,C,k3,K,e13)
    alp_24,cent_24,K24 = gauss_triple_prod_1d_xy(d4,D,l4,L,m24)
    # now we have <a1,i1,alp_13,||b2,j2,alp_24>
    # contract the (a1,i1) gaussian and (b2,j2) gaussian
    K11,alp_11,cent_11 = gauss_prod(a1,A,i1,I)
    K22,alp_22,cent_22 = gauss_prod(b2,B,j2,J)
    # this gives the alpha values and center (in 1d) of the new gaussians. Return these values and use vee in the full code
    coef1 = K13 * K11
    coef2 = K24 * K22
    return coef1, alp_11, cent_11, alp_13, cent_13, coef2, alp_22, cent_22, alp_24, cent_24

def full_contraction_3d(a1,AV,b2,BV,c3,CV,d4,DV,e13,i1,IV,j2,JV,k3,KV,l4,LV,m24):
    ndim = len(AV)
    assert ndim == len(BV)
    assert ndim == len(CV)
    assert ndim == len(DV)
    assert ndim == len(IV)
    assert ndim == len(JV)
    assert ndim == len(KV)
    assert ndim == len(LV)

    K13      = 1.0e0
    K24      = 1.0e0
    alp_11   = 0.0e0
    alp_13   = 0.0e0
    alp_22   = 0.0e0
    alp_24   = 0.0e0
    centv_11 = np.zeros(ndim)
    centv_13 = np.zeros(ndim)
    centv_22 = np.zeros(ndim)
    centv_24 = np.zeros(ndim)

    for idim in range(ndim):
        A = AV[idim]
        B = BV[idim]
        C = CV[idim]
        D = DV[idim]
        I = IV[idim]
        J = JV[idim]
        K = KV[idim]
        L = LV[idim]
        coef1, alp_11, cent_11, alp_13, cent_13, coef2, alp_22, cent_22, alp_24, cent_24 = full_contraction_1d(a1,A,b2,B,c3,C,d4,D,e13,i1,I,j2,J,k3,K,l4,L,m24)
        centv_11[idim] = cent_11
        centv_13[idim] = cent_13
        centv_22[idim] = cent_22
        centv_24[idim] = cent_24
        K13            = K13 * coef1
        K24            = K24 * coef2
    
    ans3d = K13 * K24 * vee(alp_11,centv_11,alp_13,centv_13,alp_22,centv_22,alp_24,centv_24)
    return ans3d
