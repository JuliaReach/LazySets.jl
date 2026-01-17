"""
# Extended help

    rand(::Type{HPolytope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Input

- `num_vertices` -- (optional, default: `-1`) upper bound on the number of
                    vertices of the polytope (see comment below)

### Algorithm

If `num_vertices == 0`, we create a fixed infeasible polytope.

If `num_vertices == 1`, we create a random `Singleton` and convert it.

If `dim == 1`, we create a random `Interval` and convert it.

If `dim == 2`, we create a random `VPolygon` and convert it.

Otherwise, we create a random `VPolytope` and convert it (hence also the
argument `num_vertices`). See [`rand(::Type{VPolytope})`](@ref).
"""
function rand(::Type{HPolytope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_vertices::Int=-1)
    require(@__MODULE__, :LazySets; fun_name="rand")

    rng = reseed!(rng, seed)
    if num_vertices == 0
        clist = _infeasible_constraints_list(dim; N=N)
        return HPolytope(clist)
    elseif num_vertices == 1
        P = rand(Singleton; N=N, dim=dim, rng=rng, seed=seed)
    elseif dim == 1
        if num_vertices ∉ (-1, 2)
            throw(ArgumentError("creating a 1D random polytope is only supported for 2 vertices"))
        end
        P = rand(Interval; N=N, dim=dim, rng=rng, seed=seed)
    elseif dim == 2
        P = rand(VPolygon; N=N, dim=dim, rng=rng, seed=seed, num_vertices=num_vertices)
    else
        P = rand(VPolytope; N=N, dim=dim, rng=rng, seed=seed, num_vertices=num_vertices)
    end
    return convert(HPolytope, P)
end
