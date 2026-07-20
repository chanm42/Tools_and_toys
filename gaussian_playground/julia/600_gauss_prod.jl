module GaussProd
include("700_gauss_lobe.jl")
include("800_boys_func.jl")
include("800_util_func.jl")
include("900_common_header.jl")

# Primary module for Gaussian function calculations
# Function names with 'zeroes' imply that the Gaussian centers are at zero, while more general function calculations are in the other functions

##---------------------------------------------------------##
function calc_normalization_3d(alp)
##---------------------------------------------------------##
    norm = ((2*alp) / pi)^(3/4)
    return norm
end

##---------------------------------------------------------##
function gauss_prod(a,AV,b,BV)
##---------------------------------------------------------##
    p  = (a+b)
    q  = (a*b)/(a+b)
    PV = (AV*(a/p)) + (BV*(b/p))
    sq_ab = UtilFunc.calc_dist_sq(AV,BV)
    Kab = exp(-q*sq_ab)
    return Kab,p,PV
end

##---------------------------------------------------------##
function gauss_prod_1d(alp1,AX1,beta1,BX1)
##---------------------------------------------------------##
    const_alp = (alp1*beta1) / (alp1 + beta1)
    KAB       = exp(-const_alp * (AX1 - BX1)^2)
    new_cent  = ((alp1 * AX1) + (beta1 * BX1)) / (alp1 + beta1)
    new_alp   = alp1 + beta1
    
    return KAB,new_alp,new_cent
end

##---------------------------------------------------------##
function geminal_integral_1d(alp1,AX1,gam1)
##---------------------------------------------------------##
    coef     = sqrt(pi / (alp1 + gam1))
    new_alp  = (alp1 * gam1) / (alp1 + gam1)
    new_cent = AX1
    return coef,new_alp,new_cent
end

##---------------------------------------------------------##
function gauss_trip_prod_1d(a,A,b,B,c,C)
##---------------------------------------------------------##
    gam = a + b
    P   = (a*A + b*B) / gam
    q1  = a*b / gam
    K1  = exp(-q1 * (A - B)^2)
    del = gam + c
    Q   = (gam*P + c*C) / del
    q2  = (gam*c) / del
    K2  = exp(-q2 * (P - C)^2)
    ans = K1 * K2 * sqrt(pi / (a + b + c))
    return ans
end

##---------------------------------------------------------##
function gauss_trip_prod_1d_xy(a,A,b,B,c)
##---------------------------------------------------------##
    gam1 = a + b
    q    = (a*b) / gam1
    ans  = pi * sqrt(1 / (a*c + b*c)) * exp(-q * (A - B)^2)
    return ans
end

##---------------------------------------------------------##
function gauss_trip_prod_3d_xy(a,AV,b,BV,c)
##---------------------------------------------------------##
    ndim = length(AV)
    @assert ndim == length(BV) "Dimension mismatch in Gaussian product"
    ans = 1.0e0
    for idim=1:ndim
        centA  = AV[idim]
        centB  = BV[idim]
        ans    = ans * gauss_trip_prod_1d_xy(a,centA,b,centB,c)
    end
    return ans
end

##---------------------------------------------------------##
function gauss_trip_prod_1d_xy_zeroes(a,b,c)
##---------------------------------------------------------##
    gam = a + b + c
    KAB = sqrt(pi / gam)
    new_alp = (c + c^2) / gam
    return new_alp, KAB
end

##---------------------------------------------------------##
function gauss_trip_prod_1d_xy(a,A,b,B,c)
##---------------------------------------------------------##
    gam1 = a + b
    gam2 = a + b + c
    new_alp  = (gam1 * c) / gam2
    q        = a*b / gam1
    new_cent = (a*A + b*B) / gam1
    coef     = sqrt(pi / gam2)
    KAB      = coef * exp(-q * (A-B)^2)
    return new_alp, new_cent, KAB
end

##---------------------------------------------------------##
function vee(a1,AV1,b1,BV1,c2,CV2,d2,DV2)
##---------------------------------------------------------##
    kab,pab,PVab = gauss_prod(a1,AV1,b1,BV1)
    kcd,pcd,PVcd = gauss_prod(c2,CV2,d2,DV2)
    rsq = UtilFunc.calc_dist_sq(PVab,PVcd)
    z1 = pab*pcd*rsq/(pab+pcd)
    z2 = (pab*pcd)*sqrt(pab+pcd)
    z3 = 2.0e0*sqrt(pi^5)*BoysFunc.boys_func(z1)
    ans = z3*kab*kcd/z2
    return ans
end

##---------------------------------------------------------##
function full_contraction_zeroes_1d(a1,b2,c3,d4,e13,i1,j2,k3,l4,m24)
##---------------------------------------------------------##
    alp_13,K13 = gauss_trip_prod_1d_xy_zeroes(c3,k3,e13)
    alp_24,K24 = gauss_trip_prod_1d_xy_zeroes(d4,l4,m24)
    # now we have <a1,i1,alp_13,||b2,j2,alp_24>
    # contract the (a1,i1) gaussian and (b2,j2) gaussian using 0 for the center
    K11,alp_11,cent_11 = gauss_prod(a1,0,i1,0)
    K22,alp_22,cent_22 = gauss_prod(b2,0,j2,0)
    # now calculate the vee using <alp_11,alp_13|r12|alp_22,alp_24> and multiply by the coefficients K11*K13*K22*K24
    kern       = vee(alp_11,0,alp_13,0,alp_22,0,alp_24,0)
    ans        = K13 * K11 * K24 * K22 * kern
    return ans
end

##---------------------------------------------------------##
function full_contraction_1d(a1,A,b2,B,c3,C,d4,D,e13,i1,I,j2,J,k3,K,l4,L,m24)
##---------------------------------------------------------##
    alp_13,cent_13,K13 = gauss_trip_prod_1d_xy(c3,C,k3,K,e13)
    alp_24,cent_24,K24 = gauss_trip_prod_1d_xy(d4,D,l4,L,m24)
    # now we have <a1,i1,alp_13,||b2,j2,alp_24>
    # contract the (a1,i1) gaussian and (b2,j2) gaussian
    K11,alp_11,cent_11 = gauss_prod(a1,A,i1,I)
    K22,alp_22,cent_22 = gauss_prod(b2,B,j2,J)
    # this gives the alpha values and center (in 1d) of the new gaussians. Return these values and use vee in the full code
    coef1 = K13 * K11
    coef2 = K24 * K22
    return coef1, alp_11, cent_11, alp_13, cent_13, coef2, alp_22, cent_22, alp_24, cent_24
end

##---------------------------------------------------------##
function full_contraction_3d(alp1,AV1,beta1,BV1,alp2,AV2,beta2,BV2,alp3,AV3,beta3,BV3,alp4,AV4,beta4,BV4,gam_13,gam_24)
##---------------------------------------------------------##
    ndim = length(AV1)
    @assert ndim == length(AV2)
    @assert ndim == length(AV3)
    @assert ndim == length(AV4)
    @assert ndim == length(BV1)
    @assert ndim == length(BV2)
    @assert ndim == length(BV3)
    @assert ndim == length(BV4)

    K13      = 1.0e0
    K24      = 1.0e0
    alp_11   = 0.0e0
    alp_13   = 0.0e0
    alp_22   = 0.0e0
    alp_24   = 0.0e0
    centv_11 = zeros(ndim)
    centv_13 = zeros(ndim)
    centv_22 = zeros(ndim)
    centv_24 = zeros(ndim)

    for idim=1:ndim
        AX1 = AV1[idim]
        AX2 = AV2[idim]
        AX3 = AV3[idim]
        AX4 = AV4[idim]
        BX1 = BV1[idim]
        BX2 = BV2[idim]
        BX3 = BV3[idim]
        BX4 = BV4[idim]
        
        # <g(1)g(2)g(3)g(4)g(13)|r12|h(1)h(2)h(3)h(4)h(24)>
        
        # -> <g(1)h(1)g(3)h(3)g(13)|r12|g(2)h(2)g(4)h(4)h(24)>
        
        # -> g(1)h(1) => g(11), g(3)h(3) => g(33) through gaussian product
        
        coef11, alp_11, cent_11 = gauss_prod_1d(alp1,AX1,beta1,BX1)
        coef22, alp_22, cent_22 = gauss_prod_1d(alp2,AX2,beta2,BX2)
        coef33, alp_33, cent_33 = gauss_prod_1d(alp3,AX3,beta3,BX3)
        coef44, alp_44, cent_44 = gauss_prod_1d(alp4,AX4,beta4,BX4)

        coef1 = coef11 * coef33
        coef2 = coef22 * coef44
        
        # -> <g(11)g(33)g(13)|r12|h(22)h(44)h(24)>
        
        # -> g(33)g(13) => g_3(1) through geminal integral

        coef13, alp_13, cent_13 = geminal_integral_1d(alp_33,cent_33,gam_13)
        coef24, alp_24, cent_24 = geminal_integral_1d(alp_44,cent_44,gam_24)
        
        centv_11[idim] = cent_11
        centv_13[idim] = cent_13
        centv_22[idim] = cent_22
        centv_24[idim] = cent_24

        K13            = K13 * coef1 * coef13
        K24            = K24 * coef2 * coef24
    end
            
    # -> <g(11)g_3(1)|r12|h(22)h_4(2)> => evaluate using vee
    
    ans3d = K13 * K24 * vee(alp_11,centv_11,alp_13,centv_13,alp_22,centv_22,alp_24,centv_24)
    return ans3d
end
end #module gauss_prod
