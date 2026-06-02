@validate function translate(cpa::CartesianProductArray, x::AbstractVector)
    res = Vector{LazySet}(undef, length(array(cpa)))
    s = 1
    @inbounds for (j, Xj) in enumerate(array(cpa))
        e = s + dim(Xj) - 1
        res[j] = translate(Xj, @view(x[s:e]))
        s = e + 1
    end
    return CartesianProductArray([X for X in res])
end
