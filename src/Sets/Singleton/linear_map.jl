@validate function linear_map(M::AbstractMatrix, S::Singleton)
    return Singleton(M * S.element)
end
