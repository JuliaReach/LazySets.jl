function concretize(cms::CachedMinkowskiSumArray)
    a = array(cms)
    @assert !isempty(a) "an empty Minkowski sum is not allowed"
    X = cms
    @inbounds for (i, Y) in enumerate(a)
        if i == 1
            X = concretize(Y)
        else
            X = minkowski_sum(X, concretize(Y))
        end
    end
    return X
end
