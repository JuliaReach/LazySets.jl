function extrema(P::VPolygon)
    return _extrema_vlist(P)
end

function extrema(P::VPolygon, i::Int)
    @assert 1 <= i <= dim(P) "invalid index $i for set of dimension $(dim(P))"
    return _extrema_vlist(P, i)
end
