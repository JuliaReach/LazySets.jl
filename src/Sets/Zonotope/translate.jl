function translate!(Z::Zonotope, v::AbstractVector)
    @assert length(v) == dim(Z) "cannot translate a $(dim(Z))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    Z.center .+= v
    return Z
end
