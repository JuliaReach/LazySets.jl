function translate!(E::Ellipsoid, v::AbstractVector)
    @assert length(v) == dim(E) "cannot translate a $(dim(E))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = E.center
    c .+= v
    return E
end
