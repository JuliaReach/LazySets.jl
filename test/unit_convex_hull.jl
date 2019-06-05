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

    # three vertex case in 2 dimensions
    function iscounterclockwise(result, correct_expr)
        # this function checks if the result equals any of the correct answer cyclic permutation
        result == correct_expr || result == circshift(correct_expr, 1) || result == circshift(correct_expr, 2)
    end
    ccw_points = [N[1, 1], N[-1, 1], N[-1, 0]]
    ccw_p = convex_hull!(ccw_points)
    ccw_expr = [N[1, 1], N[-1, 1], N[-1, 0]]
    @test  iscounterclockwise(ccw_p, ccw_expr) # points sorted in a counter-clockwise fashion
    cw_points = [N[-1, 1], N[1, 1], N[-1, 0]]
    cw_p = convex_hull!(cw_points)
    cw_expr = [N[1, 1], N[-1, 1], N[-1, 0]]
    @test  iscounterclockwise(cw_p, cw_expr) # points sorted in clockwise fashion
    @test ispermutation(convex_hull!([N[1, 1], N[2, 2], N[3, 3]]), [N[1, 1], N[3, 3]]) # points aligned
    @test ispermutation(convex_hull!([N[0, 1], N[0, 2], N[0, 3]]), [N[0, 1], N[0, 3]]) # three points on a vertical line
    @test convex_hull!([N[0, 1], N[0, 1], N[0, 1]]) == [N[0, 1]] # three equal points

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
