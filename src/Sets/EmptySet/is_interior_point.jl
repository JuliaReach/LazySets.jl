@validate function is_interior_point(v::AbstractVector{N}, ∅::EmptySet{N};
                                     p=N(Inf), ε=_rtol(N)) where {N<:Real}
    return false
end
