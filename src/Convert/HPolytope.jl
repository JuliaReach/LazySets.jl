# make a copy of the constraints
Base.convert(::Type{HPolytope}, P::HPolyhedron) = HPolytope(copy(constraints_list(P)))

function Base.convert(::Type{HPolytope}, ∅::EmptySet)
    clist = _infeasible_constraints_list(dim(∅); N=eltype(∅))
    return HPolytope(clist)
end
