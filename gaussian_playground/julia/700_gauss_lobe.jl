module GaussLobe
include("900_common_header.jl")
include("800_boys_func.jl")
include("800_util_func.jl")

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
function ovlap(a,AV,b,BV)
##---------------------------------------------------------##
    p   = (a+b)
    q   = (a*b)/(a+b)
    sq_ab = UtilFunc.calc_dist_sq(AV,BV)
    Kab = exp(-q*sq_ab)
    Sab = sqrt((pi/p)^3) * Kab
    return Sab
end

##---------------------------------------------------------##
function ovlap_and_ke(a,AV,b,BV)
##---------------------------------------------------------##
    p   = (a+b)
    q   = (a*b)/(a+b)
    sq_ab = UtilFunc.calc_dist_sq(AV,BV)
    Kab = exp(-q*sq_ab)
    Sab = sqrt((pi/p)^3) * Kab

    z1  = 2.0e0*q*sq_ab
    z1  = (3.0e0-z1)*q
    Tab = z1*Sab
    return Sab,Tab
end

##---------------------------------------------------------##
function vext(a,AV,b,BV,CV)
##---------------------------------------------------------##
    p   = (a+b)
    q   = (a*b)/(a+b)
    PV     = (AV*(a/p)) + (BV*(b/p))
    rsq_ab =  UtilFunc.calc_dist_sq(AV,BV)
    rsq_cp =  UtilFunc.calc_dist_sq(CV,PV)
    t1     = (2.0e0*pi/p)
    t2     = exp(-q*rsq_ab)
    t3     = p*rsq_cp
    t3     = BoysFunc.boys_func(t3)
    ans    = t1*t2*t3
    return ans
end

##---------------------------------------------------------##
function vext_soft(a,AV,b,BV,CV,r0)
##---------------------------------------------------------##
    p   = (a+b)
    q   = (a*b)/(a+b)
    PV     = (AV*(a/p)) + (BV*(b/p))
    rsq_ab =  UtilFunc.calc_dist_sq(AV,BV)
    rsq_cp =  UtilFunc.calc_dist_sq(CV,PV) + (r0^2)
    t1     = (2.0e0*pi/p)
    t2     = exp(-q*rsq_ab)
    t3     = p*rsq_cp
    t3     = BoysFunc.boys_func(t3)
    ans    = t1*t2*t3
    return ans
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
function gem_ovlap(a1,AV1,b1,BV1,c2,CV2,d2,DV2,gem)
##---------------------------------------------------------##
    kab,pab,PVab = gauss_prod(a1,AV1,b1,BV1)
    kcd,pcd,PVcd = gauss_prod(c2,CV2,d2,DV2)
    rsq = UtilFunc.calc_dist_sq(PVab,PVcd)
    z1  = (pab*pcd)+(pab*gem)+(pcd*gem)
    z2  = sqrt((pi*pi/z1)^3)
    q   = pab*pcd*gem/z1
    k12 = z2*exp(-q*rsq)
    ans = kab * kcd * k12
    return ans
end
##---------------------------------------------------------##
##---------------------------------------------------------##







end #module GaussLobe
