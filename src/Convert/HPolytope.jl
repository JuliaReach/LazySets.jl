# make a copy of the constraints
Base.convert(::Type{HPolytope}, P::HPolyhedron) = HPolytope(copy(constraints_list(P)))
