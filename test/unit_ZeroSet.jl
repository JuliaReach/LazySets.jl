Z = ZeroSet(2)
B = BallInf(ones(2), 1.)

# testing that the zero set is neutral element for the Minkowski sum
@test B ⊕ Z == B && Z ⊕ B == B

cpa = MinkowskiSumArray([B, 2.*B, 3.*B])
@test cpa ⊕ Z == cpa && Z ⊕ cpa == cpa

# test M-sum of zero set with itself
@test ZeroSet(2) ⊕ ZeroSet(2) == ZeroSet(2)

# an_element function
@test an_element(Z) ∈ Z

# subset
@test ⊆(ZeroSet(1), Singleton([0.])) && !⊆(ZeroSet(1), Singleton([2.]))
@test ⊆(ZeroSet(1), ZeroSet(1)) && !⊆(ZeroSet(1), ZeroSet(2))
