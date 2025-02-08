function cartesian_product(U1::Universe, U2::Universe)
    N = promote_type(eltype(U1), eltype(U2))
    return Universe{N}(dim(U1) + dim(U2))
end
