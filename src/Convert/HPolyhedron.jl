# make a copy of the constraints
function convert(::Type{HPolyhedron}, P::HPolytope)
    return HPolyhedron(copy(constraints_list(P)))
end

function convert(::Type{HPolyhedron}, ∅::EmptySet)
    clist = _infeasible_constraints_list(dim(∅); N=eltype(∅))
    return HPolyhedron(clist)
end
