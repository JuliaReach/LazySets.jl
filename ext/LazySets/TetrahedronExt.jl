import LazySets

using LazySets.TetrahedronModule: Tetrahedron
using LazySets.VPolytopeModule: VPolytope
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
import Base: rand
import LazySets.API: constraints_list, σ

function constraints_list(T::Tetrahedron)
    return constraints_list(convert(VPolytope, T))
end

function rand(::Type{Tetrahedron}; N::Type{<:Real}=Float64, dim::Int=3,
              rng::AbstractRNG=GLOBAL_RNG, seed::Union{Int,Nothing}=nothing)
    @assert dim == 3 "cannot create a random Tetrahedron of dimension $dim"

    rng = reseed!(rng, seed)
    P = rand(VPolytope; N=N, dim=3, rng=rng, seed=seed, num_vertices=4)
    return Tetrahedron(P.vertices)
end

"""
# Extended help

    σ(d::AbstractVector, T::Tetrahedron)

### Algorithm

This method falls back to the `VPolytope` implementation.
"""
@validate function σ(d::AbstractVector, T::Tetrahedron)
    return σ(d, convert(VPolytope, T))
end
