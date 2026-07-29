import LazySets

using LazySets: STAR, AbstractPolyhedron
using LazySets.HPolyhedronModule: HPolyhedron
using LazySets.StarModule: Star
using LazySets.UniverseModule: Universe
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: DEFAULT_COND_TOL
using ReachabilityBase.Distribution: reseed!
import Base: rand
import LazySets.API: constraints_list, isbounded, vertices_list

"""
# Extended help

    constraints_list(X::Star)

### Algorithm

See [`constraints_list(::LazySets.AbstractAffineMap)`](@ref).
"""
function constraints_list(X::Star)
    am = convert(STAR, X)
    return constraints_list(am)
end

"""
# Extended help

    isbounded(X::Star; cond_tol::Number=DEFAULT_COND_TOL)

### Algorithm

See [`isbounded(::LazySets.AbstractAffineMap)`](@ref).
"""
function isbounded(X::Star; cond_tol::Number=DEFAULT_COND_TOL)
    am = convert(STAR, X)
    return isbounded(am; cond_tol=cond_tol)
end

"""
# Extended help

    rand(::Type{Star}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Algorithm

A predicate `P` can be passed directly. If `P` is `nothing` (default), we
generate a random `HPolyhedron` of dimension `dim`.

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{Star};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              P::Union{AbstractPolyhedron,Nothing}=nothing)
    rng = reseed!(rng, seed)
    c = randn(rng, N, dim)
    if isnothing(P)
        P = rand(HPolyhedron; N=N, dim=dim, rng=rng, seed=seed)
        # may have no constraints, which is better represented as a Universe
        if isempty(P.constraints)
            P = rand(Universe; N=N, dim=dim, rng=rng, seed=seed)
        end
    end
    V = randn(rng, N, dim, LazySets.API.dim(P))
    return Star(c, V, P)
end

"""
# Extended help

    vertices_list(X::Star; apply_convex_hull::Bool=true)

### Algorithm

See [`vertices_list(::LazySets.AbstractAffineMap)`](@ref).
"""
@validate function vertices_list(X::Star; apply_convex_hull::Bool=true)
    am = convert(STAR, X)
    return vertices_list(am; apply_convex_hull=apply_convex_hull)
end
