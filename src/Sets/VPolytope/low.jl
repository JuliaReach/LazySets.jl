function low(P::VPolytope)
    return _low_vlist(P)
end

function low(P::VPolytope, i::Int)
    @assert 1 <= i <= dim(P) "invalid index $i for set of dimension $(dim(P))"
    return _low_vlist(P, i)
end