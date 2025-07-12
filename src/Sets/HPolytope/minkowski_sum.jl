@validate function minkowski_sum(P::HPolytope, Q::HPolytope;
                                 backend=nothing, algorithm=nothing, prune=true)
    res = _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
    return convert(HPolytope, res)
end
