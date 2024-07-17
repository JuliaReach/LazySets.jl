function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::Hyperplane{N},
                                 algo::AbstractLinearMapAlgorithm) where {N}
    constraints = _linear_map_hrep(M, P, algo)
    if length(constraints) == 2
        # assuming these constraints define a hyperplane
        c = first(constraints)
        return Hyperplane(c.a, c.b)
    elseif isempty(constraints)
        require(@__MODULE__, :LazySets; fun_name="linear_map")

        return Universe{N}(size(M, 1))
    else
        error("unexpected number of $(length(constraints)) constraints")
    end
end
