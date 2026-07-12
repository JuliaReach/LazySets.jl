using LazySets: default_polyhedra_backend, polyhedron
using LazySets.HPolytopeModule: HPolytope
using Polyhedra: HRep, points
import Base: convert
import LazySets.HPolytopeModule: _vertices_list

function convert(::Type{HPolytope}, P::HRep)
    return _convert_HPoly(HPolytope, P)
end

function _vertices_list(P::HPolytope; backend, prune::Bool)
    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end
    Q = polyhedron(P; backend=backend)
    if prune
        N = eltype(P)
        _removevredundancy!(Q; N=N)
    end
    return collect(points(Q))
end
