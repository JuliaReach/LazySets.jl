# convert to concrete Vector representation
function convert(::Type{HalfSpace{N,Vector{N}}},
                 hs::HalfSpace{N,<:AbstractVector{N}}) where {N}
    return HalfSpace(Vector(hs.a), hs.b)
end
