# make a copy of the constraints
function convert(::Type{HPolytope}, P::HPolyhedron)
    return HPolytope(copy(constraints_list(P)))
end

function convert(::Type{HPolytope}, ∅::EmptySet)
    clist = _infeasible_constraints_list(dim(∅); N=eltype(∅))
    return HPolytope(clist)
end
