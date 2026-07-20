module UtilFunc
include("900_common_header.jl")

const PI  = 3.141592653590e0
const N3D = 3  #number of dimensions in 3D space

function calc_dist_sq(xvec,yvec)
    ndim = length(xvec)
    distSq = 0.0e0
    for i = 1:ndim
        distSq = distSq + ((xvec[i]-yvec[i])^2)
    end
    return distSq
end

function calc_dist(xvec,yvec)
    distSq = calc_dist_sq(xvec,yvec)
    dist   = sqrt(distSq)
    return dist
end

end #module util_func