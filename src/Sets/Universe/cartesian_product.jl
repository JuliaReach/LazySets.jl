function cartesian_product(U1::Universe, U2::Universe)
    return Universe(dim(U1) + dim(U2))
end
