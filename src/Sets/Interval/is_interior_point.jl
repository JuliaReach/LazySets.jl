function is_interior_point(v::AbstractVector{N}, X::Interval{N}; p=N(Inf), ε=_rtol(N)) where {N}
    @assert ε > zero(N) "the tolerance must be strictly positive"
    @assert length(v) == 1 "a $(length(v))-dimensional vector is not compatible " *
                           "with a 1-dimensional set"
    @inbounds return min(X) + ε <= v[1] <= max(X) - ε
end
