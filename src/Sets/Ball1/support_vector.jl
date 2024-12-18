function σ(d::AbstractVector, B::Ball1)
    res = copy(B.center)
    imax = argmaxabs(d)
    res[imax] += sign(d[imax]) * B.radius
    return res
end
