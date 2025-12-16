function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::Line2D{N},
                                 algo::AbstractLinearMapAlgorithm) where {N}
    constraints = _linear_map_hrep(M, P, algo)
    if length(constraints) == 2
        # assuming these constraints define a line  # TODO assert this
        c = first(constraints)
        return Line2D(c.a, c.b)
    elseif isempty(constraints)
        return Universe{N}(size(M, 1))
    elseif length(constraints) == 4
        # projection to the origin
        S = Singleton(zeros(N, 2))
        @assert isequivalent(LazySets.HPolygon(constraints), S) "unexpected constraints"
        return S
    else
        throw(ArgumentError("unexpected number of $(length(constraints)) constraints"))
    end
end
