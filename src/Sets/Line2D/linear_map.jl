function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::Line2D{N},
                                 algo::AbstractLinearMapAlgorithm) where {N}
    constraints = _linear_map_hrep(M, P, algo)
    if length(constraints) == 2
        # assuming these constraints define a line
        c = first(constraints)
        return Line2D(c.a, c.b)
    elseif isempty(constraints)
        return Universe{N}(size(M, 1))
    else
        error("unexpected number of $(length(constraints)) constraints")
    end
end
