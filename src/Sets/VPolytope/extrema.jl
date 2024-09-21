function extrema(P::VPolytope)
    return _extrema_vlist(P)
end

function extrema(P::VPolytope, i::Int)
    @assert 1 <= i <= dim(P) "invalid index $i for set of dimension $(dim(P))"
    return _extrema_vlist(P, i)
end
