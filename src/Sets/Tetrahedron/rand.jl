function rand(::Type{Tetrahedron}; N::Type{<:Real}=Float64, rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    require(@__MODULE__, :LazySets; fun_name="rand")

    P = rand(VPolytope; N=N, dim=3, rng=rng, seed=seed, num_vertices=4)
    vertices = P.vertices
    return Tetrahedron(vertices)
end
