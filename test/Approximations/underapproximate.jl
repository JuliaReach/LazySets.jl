using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    X = Hyperrectangle(zeros(N, 2), ones(N, 2))
    U = underapproximate(X, BoxDirections{N}(2))
    @test U isa LazySet{N}
    @test U ⊆ X
end

for N in @tN([Float64, Float32])
    @static if isdefined(@__MODULE__, :Polyhedra)
        X = VPolygon([N[1, 0], N[1, 2], N[-1, 2], N[-1, 1 // 3], N[-2 // 3, 0]])
        U = underapproximate(X, Ball2)
        @test U ≈ Ball2(N[0, 1], N(1))
    end
end

for N in [Float64]
    @static if isdefined(@__MODULE__, :Ipopt)
        X = VPolygon([N[4, 2], N[2, 0], N[0, 2], N[2, 4]])
        U = underapproximate(X, Hyperrectangle)
        @test U isa Hyperrectangle{N}
        @test U ≈ Hyperrectangle(N[2, 2], ones(N, 2))
    end

    @static if isdefined(@__MODULE__, :SCS) && isdefined(@__MODULE__, :SetProg)
        @static if isdefined(@__MODULE__, :Polyhedra)
            B = Hyperrectangle(N[1, 2], N[1, 2])
            # ellipsoid from polytope (reduced error tolerance)
            E = underapproximate(B, Ellipsoid)
            E_exact = Ellipsoid(center(B), N[1 0; 0 4])
            @test E isa Ellipsoid{N}
            @test E.center ≈ E_exact.center atol = 1e-5
            @test E.shape_matrix ≈ E_exact.shape_matrix atol = 1e-4
        end

        # ellipsoid from unbounded polyhedron
        # (y position and height is rather arbitrary)
        P = HPolyhedron([HalfSpace(N[1, 0], N(2)), HalfSpace(N[-1, 0], N(0)),
                         HalfSpace(N[0, 1], N(4))])
        E = underapproximate(P, Ellipsoid; interior_point=N[1, 2])
        @test E isa Ellipsoid{N}
        @test E.center[1] ≈ N(1) atol = 1e-3
        @test ρ(N[0, 1], E) ≈ N(4) atol = 1e-3
    end
end
