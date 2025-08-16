function extrema(P::VPolygon)
    return _extrema_vlist(P)
end

@validate function extrema(P::VPolygon, i::Int)
    return _extrema_vlist(P, i)
end
