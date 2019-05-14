for N in [Float64, Rational{Int}, Float32]
    I1 = Interval(N(-1), N(1))
    I2 = Interval(N(2), N(3))
    I3 = Interval(N(-2), N(-1))
    B1 = BallInf(N[1, 1], N(2))
    B2 = BallInf(N[3, 1], N(2))
    B3 = BallInf(N[3, 3], N(2))

    # constructor
    RI1 = Rectification(I1)
    RI2 = Rectification(I2)
    RI3 = Rectification(I3)
    RB1 = Rectification(B1)
    RB2 = Rectification(B2)
    RB3 = Rectification(B3)

    # dimension
    @test dim(RI1) == dim(RI2) == dim(RI3) == 1
    @test dim(RB1) == dim(RB2) == dim(RB3) == 2
end
