import LazySets.Approximations:norm,
                               radius,
                               diameter

# approximation of a centered unit Ball2
b = Ball2([0., 0., 0.], 1.)
@test norm(b) ≈ 1.  # in the infinity norm (default)
@test radius(b) ≈ 1.  # in the infinity norm (default)
@test diameter(b) ≈ 2.  # in the infinity norm (default)

# approximations of a non-centered unit Ball2
b = Ball2([1., 2., 0.], 1.)
@test norm(b) ≈ 3.  # in the infinity norm (default)
@test radius(b) ≈ 1.  # in the infinity norm (default)
@test diameter(b) ≈ 2.  # in the infinity norm (default)
