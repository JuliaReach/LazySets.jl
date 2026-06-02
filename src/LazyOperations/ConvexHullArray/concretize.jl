function concretize(cha::ConvexHullArray)
    a = array(cha)
    if all(ispolyhedral, a)
        @assert !isempty(a) "an empty convex hull is not allowed"
        return _convex_hull_polytopes(cha)
    end

    return _concretize_lazy_array(cha)
end
