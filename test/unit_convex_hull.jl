using LazySets.Approximations: BoxDirections

for N in [Float64, Rational{Int}]


    # ====================================
    # Concrete convex hull for point sets
    # ====================================

    # corner cases in dimension 1
    @test convex_hull([Vector{N}(undef, 0)]) == [Vector{N}(undef, 0)]
    @test convex_hull([[N(0)]]) == [[N(0)]]
    @test convex_hull([[N(2)], [N(1)]]) == [[N(1)], [N(2)]]
    @test convex_hull([[N(2)], [N(2)]]) == [[N(2)]]

    # corner cases in dimension 2
    @test convex_hull([[N(0), N(0)]]) == [[N(0), N(0)]]
    @test convex_hull([[N(1), N(0)], [N(0), N(1)]]) == [[N(1), N(0)], [N(0), N(1)]]
    @test convex_hull([[N(1), N(0)], [N(1), N(0)]]) == [[N(1), N(0)]]

    # corner cases in higher dimension
    @test convex_hull([[N(0), N(0), N(0)]]) == [[N(0), N(0), N(0)]]
    # there is no sorting in higher dim
    @test convex_hull([[N(1), N(0), N(0)], [N(0), N(1), N(0)]]) == [[N(1), N(0), N(0)], [N(0), N(1), N(0)]]

    # dimension 1
    points_1D = [[N(2)], [N(1)], [N(1/2)], [N(-5)]]
    @test convex_hull(points_1D) == [[N(-5)], [N(2)]]
    convex_hull!(points_1D) # check in-place version
    @test points_1D == [[N(-5)], [N(2)]]

    # dimension 2
    points_2D = [[N(1), N(0)], [N(1), N(1)], [N(0), N(1)], [N(-1), N(-1)], [N(1/2), N(1/2)]]
    @test convex_hull(points_2D) == [[N(-1), N(-1)], [N(1), N(0)], [N(1), N(1)], [N(0), N(1)]]
    convex_hull!(points_2D) # check in-place version
    @test points_2D == [[N(-1), N(-1)], [N(1), N(0)], [N(1), N(1)], [N(0), N(1)]]

    # higher dimension
    if test_suite_polyhedra && N != Float32 # TODO : check why it breaks with Float32
        points_3D = [[N(1), N(0), N(4)], [N(1), N(1), N(5)], [N(0), N(1), N(6)],
                     [N(-1), N(-1), N(-7)], [N(1/2), N(1/2), N(-8)], [N(0), N(0), N(0)],
                     [N(1), N(2), N(3)]]
        @test convex_hull(points_3D) == ispermutation([[N(1), N(0), N(4)], [N(1), N(1), N(5)],
                                         [N(0), N(1), N(6)], [N(-1), N(-1), N(-7)],
                                         [N(1/2), N(1/2), N(-8)], [N(1), N(2), N(3)]])
        convex_hull!(points_3D) # check in-place version
        @test points_3D == ispermutation([[N(1), N(0), N(4)], [N(1), N(1), N(5)],
                            [N(0), N(1), N(6)], [N(-1), N(-1), N(-7)],
                            [N(1/2), N(1/2), N(-8)], [N(1), N(2), N(3)]])
    end

    # ============================
    # Binary concrete convex hull 
    # ============================
    V1 = VPolygon([[N(1), N(0)], [N(1), N(1)], [N(0), N(1)]])
    V2 = VPolygon([[N(-1), N(-1)], [N(1/2), N(1/2)]])
    ch = convex_hull(V1, V2)
    @test vertices_list(ch) == ispermutation([[N(-1), N(-1)], [N(1), N(0)], [N(1), N(1)], [N(0), N(1)]])

    if test_suite_polyhedra && N != Float32 # TODO : check why it breaks with Float32
        V1 = VPolytope([[N(1), N(0), N(4)], [N(1), N(1), N(5)], [N(0), N(1), N(6)]])
        V2 = VPolytope([[N(-1), N(-1), N(-7)], [N(1/2), N(1/2), N(-8)], [N(0), N(0), N(0)],
                       [N(1), N(2), N(3)]])
        ch = convex_hull(V1, V2)
        @test ch isa VPolytope
        @test vertices_list(ch) == ispermutation([[N(1), N(0), N(4)], [N(1), N(1), N(5)],
                                    [N(0), N(1), N(6)], [N(-1), N(-1), N(-7)],
                                    [N(1/2), N(1/2), N(-8)], [N(1), N(2), N(3)]])
    end
end
