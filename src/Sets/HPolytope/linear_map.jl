function _linear_map_hrep_helper(M::AbstractMatrix, P::HPolytope,
                                 algo::AbstractLinearMapAlgorithm)
    constraints = _linear_map_hrep(M, P, algo)
    return HPolytope(constraints)
end
