import LazySets.Approximations.ballinf_approximation

for N in [Float64, Float32, Rational{Int}]
    # BallInf approximation of a 3D unit ball in the 1-norm centered at [1,2,0]
    b = Ball1(N[1, 2, 0], N(1))
    bi = ballinf_approximation(b)
    biexp = BallInf(N[1, 2, 0], N(1))
    @test bi.center ≈ biexp.center
    @test bi.radius ≈ biexp.radius

    # BallInf approximation of a 2D polygon whose vertices are
    # (0,1), (1,1), (2,-1), (-1,0)
    p = HPolygon{N}()
    addconstraint!(p, LinearConstraint(N[0, 1], N(1)))
    addconstraint!(p, LinearConstraint(N[-1, 1], N(1)))
    addconstraint!(p, LinearConstraint(N[-1, -3], N(1)))
    addconstraint!(p, LinearConstraint(N[2, 1], N(3)))
    bi = ballinf_approximation(p)
    biexp = BallInf(N[0.5, 0], N(1.5))
    @test bi.center ≈ biexp.center
    @test bi.radius ≈ biexp.radius
end
