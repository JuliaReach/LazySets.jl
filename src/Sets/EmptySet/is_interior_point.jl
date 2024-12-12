function is_interior_point(v::AbstractVector{N}, ∅::EmptySet{N}; p=N(Inf),
                           ε=_rtol(N)) where {N<:Real}
    @assert length(v) == dim(∅) "incompatible dimensions $(length(v)) and $(dim(∅))"
    @assert ε > zero(N) "the tolerance must be strictly positive"

    return false
end
