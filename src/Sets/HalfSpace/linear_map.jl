function _linear_map_hrep_helper(M::AbstractMatrix, hs::HalfSpace,
                                 algo::AbstractLinearMapAlgorithm)
    constraints = _linear_map_hrep(M, hs, algo)
    if length(constraints) == 1
        return first(constraints)
    elseif isempty(constraints)
        require(@__MODULE__, :LazySets; fun_name="linear_map")

        N = promote_type(eltype(M), eltype(hs))
        return Universe{N}(size(M, 1))
    else
        require(@__MODULE__, :LazySets; fun_name="linear_map")

        return HPolyhedron(constraints)
    end
end
