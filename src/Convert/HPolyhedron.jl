# make a copy of the constraints
Base.convert(::Type{HPolyhedron}, P::HPolytope) = HPolyhedron(copy(constraints_list(P)))
