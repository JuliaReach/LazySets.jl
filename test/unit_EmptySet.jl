E = EmptySet()
B = BallInf(ones(2), 1.)

# testing that the empty set is an absorbing element for the cartesian product
@test B * E isa EmptySet && E * B isa EmptySet
# testing the mathematical alias ×
@test B × E isa EmptySet && E × B isa EmptySet
# test ∅ alias
@test B × E isa ∅

cpa = CartesianProductArray([B, 2.*B, 3.*B])
@test cpa * E isa EmptySet && E * cpa isa EmptySet
@test cpa × E isa EmptySet && E × cpa isa EmptySet

# testing that the empty set is an absorbing element for the Minkowski sum
@test B + E isa EmptySet && B + E isa EmptySet
# testing the mathematical alias ⊕ 
@test B ⊕ E isa EmptySet && B ⊕ E isa EmptySet

msa = MinkowskiSumArray([B, 2.*B, 3.*B])
@test msa + E isa EmptySet && E + msa isa EmptySet
@test msa ⊕ E isa EmptySet && E ⊕ msa isa EmptySet
