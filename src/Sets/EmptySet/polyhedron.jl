function polyhedron(∅::EmptySet; backend=default_polyhedra_backend(∅))
    N = eltype(∅)
    P = Polyhedra.vrep(Vector{N}[]; d=dim(∅))
    return Polyhedra.polyhedron(P, backend)
end
