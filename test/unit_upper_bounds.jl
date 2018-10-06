import LazySets.Approximations.ρ_upper_bound

X = Hyperrectangle(low=[0.0, 0.0], high=[1.0, 1.0])
d = [1.0, 2.0]
@test ρ(d, X) == ρ_upper_bound(d, X)

using Optim
@test ρ(d, X) == ρ_upper_bound(d, X ∩ X)
M = [1.0 1.0; 1.0 1.0]
@test isapprox(ρ_upper_bound(d, X ∩ (M * X)), 3.0, atol=1e-6)
