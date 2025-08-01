using Test, LazySets
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
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

    # empty set
    E = EmptySet{N}(2)
    @test ballinf_approximation(E) == E

    if N == Float64
        # robustness (see issue #2532): the set has two contradicting constraints
        # => the set is empty, but requires to use high precision
        # x >= 1.0000000000000002
        # x <= 1
        s1 = HalfSpace([-1.0], -1.0000000000000002)
        s2 = HalfSpace([1.0], 1.0)
        P = s1 ∩ s2
        @test ballinf_approximation(P) == BallInf([1.0], 0.0)

        s1big = HalfSpace([-big(1.0)], -big(1.0000000000000002))
        s2big = HalfSpace([big(1.0)], big(1.0))
        Pbig = HPolytope([s1big, s2big])
        @test ballinf_approximation(Pbig) == EmptySet{BigFloat}(1)
    end
end
