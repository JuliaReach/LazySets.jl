function rand(::Type{Tetrahedron}; N::Type{<:Real}=Float64, dim::Int=3,
              rng::AbstractRNG=GLOBAL_RNG, seed::Union{Int,Nothing}=nothing)
    @assert dim == 3 "cannot create a random Tetrahedron of dimension $dim"
    require(@__MODULE__, :LazySets; fun_name="rand")

    rng = reseed!(rng, seed)
    P = rand(VPolytope; N=N, dim=3, rng=rng, seed=seed, num_vertices=4)
    return Tetrahedron(P.vertices)
end
