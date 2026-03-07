@validate function exponential_map(M::AbstractMatrix, U::Universe)
    return _linear_map_universe_invertible(M, U)
end
