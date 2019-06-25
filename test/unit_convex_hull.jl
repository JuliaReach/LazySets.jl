for N in [Float64, Rational{Int}]
    # ====================================
    # Concrete convex hull for point sets
    # ====================================

    # corner cases in dimension 1
    @test convex_hull([Vector{N}(undef, 0)]) == [Vector{N}(undef, 0)]
    @test convex_hull([[N(0)]]) == [[N(0)]]
    @test convex_hull([[N(2)], [N(1)]]) == [[N(2)], [N(1)]]
    @test convex_hull([[N(2)], [N(2)]]) == [[N(2)]]

    # corner cases in dimension 2
    @test convex_hull([[N(0), N(0)]]) == [[N(0), N(0)]]
    @test convex_hull([[N(1), N(0)], [N(0), N(1)]]) == [[N(1), N(0)], [N(0), N(1)]]
    @test convex_hull([[N(1), N(0)], [N(1), N(0)]]) == [[N(1), N(0)]]

    # test corner cases with one and two vectors (see #876)
    p1 = [1., 2.]
    p2 = [1., 3.]
    @test convex_hull([p1]) == [p1]
    @test convex_hull([p1, p2]) == [p1, p2]
    @test convex_hull([p2, p1]) == [p2, p1] # no sorting

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

    # three vertex case in 2 dimensions
    ccw_points = [N[1, 1], N[-1, 1], N[-1, 0]]
    ccw_p = convex_hull!(ccw_points)
    ccw_expr = [N[1, 1], N[-1, 1], N[-1, 0]]
    @test is_cyclic_permutation(ccw_p, ccw_expr) # points sorted in a counter-clockwise fashion
    cw_points = [N[-1, 1], N[1, 1], N[-1, 0]]
    cw_p = convex_hull!(cw_points)
    cw_expr = [N[1, 1], N[-1, 1], N[-1, 0]]
    @test is_cyclic_permutation(cw_p, cw_expr) # points sorted in clockwise fashion
    @test ispermutation(convex_hull!([N[1, 1], N[2, 2], N[3, 3]]), [N[1, 1], N[3, 3]]) # points aligned
    @test ispermutation(convex_hull!([N[0, 1], N[0, 2], N[0, 3]]), [N[0, 1], N[0, 3]]) # three points on a vertical line
    @test convex_hull!([N[0, 1], N[0, 1], N[0, 1]]) == [N[0, 1]] # three equal points

    # four vertex case in 2 dimentions
    A = N[1, 0]
    B = N[1, 1]
    C = N[-1, 1]
    D = N[-1, 0]
    expr = [A, B, C, D]
    points = [A, B, C, D]
    @test is_cyclic_permutation(convex_hull!(points), expr) # ABCD
    points = [A, D, C, B]
    @test is_cyclic_permutation(convex_hull!(points), expr) # ADCB
    points = [A, B, D, C]
    @test is_cyclic_permutation(convex_hull!(points), expr) # ABDC
    points = [A, D, B, C]
    @test is_cyclic_permutation(convex_hull!(points), expr) # ADBC
    points = [D, A, C, B]
    @test is_cyclic_permutation(convex_hull!(points), expr) # DACB
    points = [A, C, D, B]
    @test is_cyclic_permutation(convex_hull!(points), expr) # ACDB
    A = N[0, 1]
    B = N[-1, -1]
    C = N[1, -1]
    D = N[0, 0]
    expr = [A, B, C]
    points = [A, B, C, D]
    @test is_cyclic_permutation(convex_hull!(points), expr) # ABC
    points = [A, C, B, D]
    @test is_cyclic_permutation(convex_hull!(points), expr) # CBA
    points = [D, B, C, A]
    @test is_cyclic_permutation(convex_hull!(points), expr) # BCD
    points = [D, C, B, A]
    @test is_cyclic_permutation(convex_hull!(points), expr) # DCB
    points = [B, D, C, A]
    @test is_cyclic_permutation(convex_hull!(points), expr) # ACD
    points = [C, D, B, A]
    @test is_cyclic_permutation(convex_hull!(points), expr) # CAD
    points = [B, C, D, A]
    @test is_cyclic_permutation(convex_hull!(points), expr) # ABD
    points = [C, B, D, A]
    @test is_cyclic_permutation(convex_hull!(points), expr) # ADB
    @test ispermutation(convex_hull!([N[1, 1], N[2, 2], N[3, 3], N[4, 4]]), [N[1, 1], N[4, 4]]) # points aligned

    # five-vertices case in 2D
    points = to_N(N, [[0.9, 0.2], [0.4, 0.6], [0.2, 0.1], [0.1, 0.3], [0.3, 0.28]])
    points_copy = copy(points)
    sorted = to_N(N, [[0.1, 0.3], [0.2, 0.1], [0.9, 0.2], [0.4, 0.6]])
    @test ispermutation(convex_hull!(points, algorithm="monotone_chain"), sorted)
    @test ispermutation(convex_hull!(points, algorithm="monotone_chain_sorted"), sorted)
    @test_throws ErrorException convex_hull!(points_copy, algorithm="")

    # higher dimension
    if test_suite_polyhedra && N != Float32 # no backend supporting Float32
        points_3D = [[N(1), N(0), N(4)], [N(1), N(1), N(5)], [N(0), N(1), N(6)],
                     [N(-1), N(-1), N(-7)], [N(1/2), N(1/2), N(-8)], [N(0), N(0), N(0)],
                     [N(1), N(2), N(3)]]
        @test ispermutation(convex_hull(points_3D), [[N(1), N(0), N(4)], [N(1), N(1), N(5)],
                                         [N(0), N(1), N(6)], [N(-1), N(-1), N(-7)],
                                         [N(1/2), N(1/2), N(-8)], [N(1), N(2), N(3)]])
        convex_hull!(points_3D) # check in-place version
        @test ispermutation(points_3D, [[N(1), N(0), N(4)], [N(1), N(1), N(5)],
                            [N(0), N(1), N(6)], [N(-1), N(-1), N(-7)],
                            [N(1/2), N(1/2), N(-8)], [N(1), N(2), N(3)]])
    end

    # ============================
    # Binary concrete convex hull
    # ============================

    V1 = VPolygon([[N(1), N(0)], [N(1), N(1)], [N(0), N(1)]])
    V2 = VPolygon([[N(-1), N(-1)], [N(1/2), N(1/2)]])
    ch = convex_hull(V1, V2)
    @test ispermutation(vertices_list(ch), [[N(-1), N(-1)], [N(1), N(0)], [N(1), N(1)], [N(0), N(1)]])

    if test_suite_polyhedra && N != Float32 # no backend supporting Float32
        V1 = VPolytope([[N(1), N(0), N(4)], [N(1), N(1), N(5)], [N(0), N(1), N(6)]])
        V2 = VPolytope([[N(-1), N(-1), N(-7)], [N(1/2), N(1/2), N(-8)], [N(0), N(0), N(0)],
                       [N(1), N(2), N(3)]])
        ch = convex_hull(V1, V2)
        @test ch isa VPolytope
        @test ispermutation(vertices_list(ch), [[N(1), N(0), N(4)], [N(1), N(1), N(5)],
                                    [N(0), N(1), N(6)], [N(-1), N(-1), N(-7)],
                                    [N(1/2), N(1/2), N(-8)], [N(1), N(2), N(3)]])
    end
end
