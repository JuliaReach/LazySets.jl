@validate function exponential_map(M::AbstractMatrix, X::Interval)
    @inbounds e = exp(M[1, 1])
    return Interval(e * X.dat)
end
