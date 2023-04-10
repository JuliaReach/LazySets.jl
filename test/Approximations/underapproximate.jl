for N in [Float64, Rational{Int}, Float32]
    X = Hyperrectangle(zeros(N, 2), ones(N, 2))
    U = underapproximate(X, BoxDirections{N}(2))
    @test U ⊆ X

    X = VPolygon([N[4, 2], N[2, 0], N[0, 2], N[2, 4]])
    U = underapproximate(X, Hyperrectangle)
    @test U ≈ Hyperrectangle(N[2, 2], ones(N, 2))
end

for N in [Float64, Float32]
    X = VPolygon([N[1, 0], N[1, 2], N[-1, 2], N[-1, 1//3], N[-2//3, 0]])
    U = underapproximate(X, Ball2)
    @test U ≈ Ball2(N[0, 1], N(1))
end

for N in [Float64]
    # ellipsoid from polytope (reduced error tolerance)
    B = Hyperrectangle(N[1, 2], N[1, 2])
    E = underapproximate(B, Ellipsoid)
    E_exact = Ellipsoid(center(B), N[1 0; 0 4])
    @test E isa Ellipsoid{N}
    @test E.center ≈ E_exact.center atol=1e-5
    @test E.shape_matrix ≈ E_exact.shape_matrix atol=1e-4
    # ellipsoid from unbounded polyhedron
    # (y position and height is rather arbitrary)
    P = HPolyhedron([HalfSpace(N[1, 0], N(2)), HalfSpace(N[-1, 0], N(0)),
                     HalfSpace(N[0, 1], N(4))])
    E = underapproximate(P, Ellipsoid; interior_point=center(B))
    @test E isa Ellipsoid{N}
    @test E.center[1] ≈ N(1) atol=1e-3
    @test ρ(N[0, 1], E) ≈ N(4) atol=1e-3
end
