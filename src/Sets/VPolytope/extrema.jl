function extrema(P::VPolytope)
    return _extrema_vlist(P)
end

@validate function extrema(P::VPolytope, i::Int)
    return _extrema_vlist(P, i)
end
