function polyhedron(∅::EmptySet; backend=default_polyhedra_backend(∅))
    N = eltype(∅)
    P = Polyhedra.vrep(Vector{N}[]; d=dim(∅))
    return polyhedron(P, backend)
end
