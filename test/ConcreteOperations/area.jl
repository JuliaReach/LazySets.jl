for N in @tN([Float64, Float32, Rational{Int}])
    # not implemented for dimension other than 2 or 3
    for d in (1, 4)
        p = BallInf(zeros(N, d), N(1))
        @test_throws DimensionMismatch area(p)
    end

    # sets with zero area
    p = Singleton(N[0, 1])
    @test area(p) ≈ N(0)
    p = Interval(N(0), N(1)) × Interval(N(0), N(0))
    @test area(p) ≈ N(0)

    # triangle
    p = VPolygon([N[13, 14], N[16, 30], N[50, 10]])
    @test area(p) ≈ N(302)

    # quadrilateral
    vlist = [N[3, 4], N[5, 11], N[12, 8], N[9, 5], N[5, 6]]
    p = VPolygon(vlist) # convex hull => 4 points
    @test area(p) ≈ N(35)

    # more than four vertices
    p2 = VPolygon([N[1, 0], N[0, 0], N[0, 1], N[1 / 2, 2], [1, 1]])
    @test area(p2) == N(3 / 2)

    # non-convex polygon (via list of vertices)
    @test LazySets._area_polygon(vlist) ≈ N(30)
end
