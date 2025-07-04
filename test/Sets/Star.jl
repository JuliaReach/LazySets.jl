using LazySets, Test
using LazySets.ReachabilityBase.Arrays: ispermutation

for N in [Float64, Float32, Rational{Int}]
    # constructor with basis matrix
    S = Star(N[3, 3], N[1 0; 0 1], BallInf(N[0, 0], N(1)))

    # constructor vector-of-vectors basis
    W = Star(N[3, 3], [N[1, 0], N[0, 1]], BallInf(N[0, 0], N(1)))

    @test S ≈ W

    # getter functions
    @test center(S) == N[3, 3]
    @test center(S, 1) == center(S, 2) == N(3)
    @test basis(S) == N[1 0; 0 1]
    @test predicate(S) == BallInf(N[0, 0], N(1))
    @test dim(S) == 2

    # support function
    @test ρ(N[1, 0], S) == N(4)
    @test -ρ(N[-1, 0], S) == N(2)
    @test ρ(N[0, 1], S) == N(4)
    @test -ρ(N[0, -1], S) == N(2)

    # support vector
    @test σ(N[1, 1], S) == N[4, 4]
    @test σ(N[-1, -1], S) == N[2, 2]
    @test σ(N[1, -1], S) == N[4, 2]
    @test σ(N[-1, 1], S) == N[2, 4]

    # conversion from a polyhedral set
    P = HPolygon([HalfSpace(N[1, 1], N(-0.5)),
                  HalfSpace(N[-1.5, 0.5], N(4)),
                  HalfSpace(N[-0.5, -0.05], N(1.5)),
                  HalfSpace(N[1.5, -0.25], N(-1.5))])
    S = convert(Star, P)
    @test basis(S) == N[1 0; 0 1]
    @test center(S) == N[0, 0]
    @test predicate(S) == P

    # ispolyhedral
    @test ispolyhedral(S)

    # intersection with HalfSpace
    H = HalfSpace(N[0, 1], N(0))
    I = intersection(S, H)
    intersection!(S, H)
    addconstraint!(P, H)
    @test isequivalent(I, P)
    @test isequivalent(S, P)

    # check that we can intersect polyhedra that are axis-aligned
    B = BallInf(N[0, 0], N(1))
    S = Star(N[0, 0], N[1 0; 0 1], B)
    I = intersection(S, H)
    @test isequivalent(I, intersection(B, H))

    # concrete linear map (returns a star)
    M = N[2 0; 0 2]
    lm = linear_map(M, S)
    @test isa(lm, Star)
    @test isequivalent(lm, Star(N[0, 0], M, B))

    # concrete affine map (returns a star)
    v = N[1, 2]
    am = affine_map(M, S, v)
    @test isa(am, Star)
    @test isequivalent(am, Star(v, M, B))

    # an_element and ∈
    e = an_element(S)
    @test e ∈ B && e ∈ S

    # isempty
    @test !isempty(S)

    # isbounded
    @test isbounded(S)

    # vertices_list
    @test ispermutation(vertices_list(S), vertices_list(B))
end

for N in [Float64]
    # random star
    rand(Star)
end

# isoperationtype
@test !isoperationtype(Star)
