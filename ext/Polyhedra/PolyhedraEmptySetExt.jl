using LazySets: default_polyhedra_backend, dim
using LazySets.EmptySetModule: EmptySet
using Polyhedra: vrep
import LazySets: polyhedron

function polyhedron(∅::EmptySet; backend=default_polyhedra_backend(∅))
    N = eltype(∅)
    P = vrep(Vector{N}[]; d=dim(∅))
    return Polyhedra.polyhedron(P, backend)
end
