function load_Polyhedra_polyhedron()
    return quote
        import .Polyhedra: polyhedron

        """
            polyhedron(P::VPolytope;
                       [backend]=default_polyhedra_backend(P),
                       [relative_dimension]=nothing)

        Return a `VRep` polyhedron from `Polyhedra.jl` given a polytope in vertex
        representation.

        ### Input

        - `P`       -- polytope in vertex representation
        - `backend` -- (optional, default: `default_polyhedra_backend(P)`) the
                       backend for polyhedral computations; see [Polyhedra's
                       documentation](https://juliapolyhedra.github.io/) for further
                       information
        - `relative_dimension` -- (default, optional: `nothing`) an integer representing
                                  the (relative) dimension of the polytope; this
                                  argument is mandatory if the polytope is empty

        ### Output

        A `VRep` polyhedron.

        ### Notes

        The *relative dimension* (or just *dimension*) refers to the dimension of the
        set relative to itself, independently of the ambient dimension. For example, a
        point has (relative) dimension zero, and a line segment has (relative) dimension
        one.

        In this library, `dim` always returns the ambient dimension of the set,
        such that a line segment in two dimensions has dimension two. However,
        `Polyhedra.dim` will assign a dimension equal to one to a line segment
        because it uses a different convention.
        """
        function polyhedron(P::VPolytope;
                            backend=default_polyhedra_backend(P),
                            relative_dimension=nothing)
            if isempty(P)
                if isnothing(relative_dimension)
                    throw(ArgumentError("the conversion to a `Polyhedra.polyhedron` requires the " *
                                        "(relative) dimension of the `VPolytope` to be known, " *
                                        "but it cannot be inferred from an empty set; use the " *
                                        "keyword argument `relative_dimension`"))
                end
                return polyhedron(Polyhedra.vrep(P.vertices; d=relative_dimension), backend)
            end
            return polyhedron(Polyhedra.vrep(P.vertices), backend)
        end
    end
end  # load_Polyhedra_polyhedron
