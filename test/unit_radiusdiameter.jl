for N in [Float64, Rational{Int}, Float32]
    # =======
    #  Ball1
    # =======

    # approximation of a centered unit Ball1
    b = Ball1(N[0., 0., 0.], N(1.))
    @test norm(b) ≈ N(1.)  # in the infinity norm (default)
    @test radius(b) ≈ N(1.)  # in the infinity norm (default)
    @test diameter(b) ≈ N(2.)  # in the infinity norm (default)

    # approximations of a non-centered unit Ball1
    b = Ball1(N[1., 2., 0.], N(1.))
    @test norm(b) ≈ N(3.)  # in the infinity norm (default)
    @test radius(b) ≈ N(1.)  # in the infinity norm (default)
    @test diameter(b) ≈ N(2.)  # in the infinity norm (default)


    # ================
    #  Hyperrectangle
    # ================

    # metrics in the infinity norm (default)
    h = Hyperrectangle(N[-1., 2.], N[0.2, 0.5])
    @test norm(h, Inf) ≈ N(2.5)
    @test radius(h, Inf) ≈ N(0.5)
    @test diameter(h, Inf) ≈ N(1.)

    # metrics in the 2-norm
    @test norm(h, 2) ≈ N(sqrt(2.5^2 + 1.2^2))
    @test radius(h, 2) ≈ N(sqrt(.5^2 + .2^2))
    @test diameter(h, 2) ≈ N(2*sqrt(.5^2 + .2^2))

    # =========
    #  BallInf
    # =========

    # metrics in the infinity norm (default)
    b = BallInf(N[-1., -2., -4], N(0.2))
    @test norm(b, Inf) ≈ N(4.2)
    @test radius(b, Inf) ≈ N(0.2)
    @test diameter(b, Inf) ≈ N(0.4)

    # metrics in the 2-norm
    @test norm(b, 2) ≈ N(sqrt(1.2^2 + 2.2^2 + 4.2^2))
    @test radius(b, 2) ≈ N(sqrt(0.2^2 * 3))
    @test diameter(b, 2) ≈ N(2*sqrt(0.2^2 * 3))
end
