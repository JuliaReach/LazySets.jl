function exponential_map(M::AbstractMatrix, X::Interval)
    @assert size(M) == (1, 1) "cannot apply an exponential map of dimension " *
                              "$(size(M)) to an interval"
    @inbounds e = exp(M[1, 1])
    return Interval(e * X.dat)
end
