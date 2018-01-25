import LazySets.Approximations.overapproximate

for N in [Float64, Float32] # TODO Rational{Int}
    # Approximation of a 2D centered unit ball in norm 1
    # All vertices v should be like this:
    # ‖v‖ >= 1 and ‖v‖ <= 1+ɛ
    # Where ɛ is the given error bound
    b = Ball1(N[0., 0.], N(1.))
    ɛ = N(.01)
    p = tovrep(overapproximate(b, ɛ))
    for v in p.vertices_list
    @test norm(v) >= N(1.)
    @test norm(v) <= N(1.+ɛ)
    end

    # Check that there are no redundant constraints for a ballinf
    b = BallInf(N[0.5, 0.5], N(0.1))
    lcl = overapproximate(b, N(.001)).constraints_list
    @test length(lcl) == 4
    @test lcl[1].a == N[1.0, 0.0]
    @test lcl[1].b == N(0.6)
    @test lcl[2].a == N[0.0, 1.0]
    @test lcl[2].b == N(0.6)
    @test lcl[3].a == N[-1.0, 0.0]
    @test lcl[3].b == N(-0.4)
    @test lcl[4].a == N[0.0, -1.0]
    @test lcl[4].b == N(-0.4)

    # Check that there are no redundant constraints for a HPolygon (octagon)
    p = HPolygon{N}()
    addconstraint!(p, LinearConstraint(N[1.0, 0.0], N(1.0)))
    addconstraint!(p, LinearConstraint(N[0.0, 1.0], N(1.0)))
    addconstraint!(p, LinearConstraint(N[-1.0, 0.0], N(1.0)))
    addconstraint!(p, LinearConstraint(N[0.0, -1.0], N(1.0)))
    addconstraint!(p, LinearConstraint(N[sqrt(2.0)/2.0, sqrt(2.0)/2.0], N(1.0)))
    addconstraint!(p, LinearConstraint(N[-sqrt(2.0)/2.0, sqrt(2.0)/2.0], N(1.0)))
    addconstraint!(p, LinearConstraint(N[sqrt(2.0)/2.0, -sqrt(2.0)/2.0], N(1.0)))
    addconstraint!(p, LinearConstraint(N[-sqrt(2.0)/2.0, -sqrt(2.0)/2.0], N(1.0)))
    lcl = overapproximate(p, N(.001)).constraints_list
    @test length(lcl) == 8
    @test lcl[1].a ≈ N[1.0, 0.0]
    @test lcl[1].b ≈ N(1.0)
    @test lcl[2].a ≈ N[sqrt(2.0)/2.0, sqrt(2.0)/2.0]
    @test lcl[2].b ≈ N(1.0)
    @test lcl[3].a ≈ N[0.0, 1.0]
    @test lcl[3].b ≈ N(1.0)
    @test lcl[4].a ≈ N[-sqrt(2.0)/2.0, sqrt(2.0)/2.0]
    @test lcl[4].b ≈ N(1.0)
    @test lcl[5].a ≈ N[-1.0, 0.0]
    @test lcl[5].b ≈ N(1.0)
    @test lcl[6].a ≈ N[-sqrt(2.0)/2.0, -sqrt(2.0)/2.0]
    @test lcl[6].b ≈ N(1.0)
    @test lcl[7].a ≈ N[0.0, -1.0]
    @test lcl[7].b ≈ N(1.0)
    @test lcl[8].a ≈ N[sqrt(2.0)/2.0, -sqrt(2.0)/2.0]
    @test lcl[8].b ≈ N(1.0)
end
