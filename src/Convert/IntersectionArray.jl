"""
    convert(::Type{IntersectionArray}, P::AbstractPolytope{N}) where {N<:AbstractFloat}

Represent a nonempty, bounded, positive-dimensional polytope as an exact
intersection of zonotopes. Each nonzero constraint normal defines one
containing zonotope whose supporting slab enforces that constraint, up to
the numerical precision of the polytope's support function.
"""
function convert(::Type{IntersectionArray}, P::AbstractPolytope{N}) where {N<:AbstractFloat}
    n = dim(P)
    if n < 1 || isempty(P) || !isbounded(P)
        throw(ArgumentError("expected a nonempty bounded polytope of positive dimension"))
    end

    # Reuse coordinate bounds for every slab; only its normal needs new bounds.
    lower, upper = extrema(P)
    zonotopes = Zonotope{N,Vector{N},Matrix{N}}[]
    for H in constraints_list(P)
        a = H.a
        all(iszero, a) && continue
        pivot = argmax(abs.(a))
        transform = Matrix{N}(I, n, n)
        transform[pivot, :] = a
        lo = Vector{N}(lower)
        hi = Vector{N}(upper)
        lo[pivot] = -ρ(-a, P)
        hi[pivot] = H.b
        midpoint = lo / 2 + hi / 2
        radii = hi / 2 - lo / 2
        push!(zonotopes, Zonotope(transform \ midpoint,
                                  transform \ Matrix(Diagonal(radii))))
    end
    return IntersectionArray(zonotopes)
end
