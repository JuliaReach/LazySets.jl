# make a copy of the constraints
function convert(::Type{HPolyhedron}, P::HPolytope)
    return HPolyhedron(copy(constraints_list(P)))
end
