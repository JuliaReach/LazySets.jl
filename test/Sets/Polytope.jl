using LazySets, Test, LinearAlgebra
using LazySets: linear_map_inverse, affine_map_inverse
using LazySets.ReachabilityBase.Arrays: SingleEntryVector, ispermutation

for N in [Float64, Float32, Rational{Int}]
    # random polytopes
    @test_throws ArgumentError rand(HPolytope; N=N, dim=1, num_vertices=3)
    p = rand(HPolytope; N=N, dim=1, num_vertices=0)
    @test p isa HPolytope{N} && dim(p) == 1 && isempty(p)
    p = rand(VPolytope; N=N, num_vertices=0)
    @test p isa VPolytope{N} && length(p.vertices) == 0

    # -----
    # H-rep
    # -----

    # constructor from matrix and vector
    A = [N(1) N(2); N(-1) N(1)]
    b = [N(1), N(2)]
    p = HPolytope(A, b)
    c = p.constraints
    @test c isa Vector{LinearConstraint{N,Vector{N}}}
    @test c[1].a == N[1, 2] && c[1].b == N(1)
    @test c[2].a == N[-1, 1] && c[2].b == N(2)
    @test_throws AssertionError HPolytope(N[1 0; 0 1], N[1, 1];
                                          check_boundedness=true)

    # convert back to matrix and vector
    A2, b2 = tosimplehrep(p)
    @test A == A2 && b == b2

    # 2D polytope
    p = HPolytope{N}()
    c1 = LinearConstraint(N[2, 2], N(12))
    c2 = LinearConstraint(N[-3, 3], N(6))
    c3 = LinearConstraint(N[-1, -1], N(0))
    c4 = LinearConstraint(N[2, -4], N(0))
    addconstraint!(p, c3)
    addconstraint!(p, c1)
    addconstraint!(p, c4)
    addconstraint!(p, c2)

    # support vector of polytope with no constraints
    @test_throws ErrorException σ(N[0], HPolytope{N}())

    # boundedness
    @test isbounded(p) && isbounded(p, false) && isboundedtype(typeof(p))
    p2 = HPolytope{N}()
    @test isbounded(p2) && !isbounded(p2, false) && isboundedtype(typeof(p2))

    # ispolyhedral
    @test ispolyhedral(p)

    # isuniversal
    answer, w = isuniversal(p, true)
    @test !isuniversal(p) && !answer && w ∉ p

    # membership
    @test N[5 / 4, 7 / 4] ∈ p
    @test N[4, 1] ∉ p

    # singleton list
    @test length(singleton_list(p)) == 4

    @static if isdefined(@__MODULE__, :Polyhedra)
        # conversion to and from Polyhedra's VRep data structure
        cl = constraints_list(convert(HPolytope, polyhedron(p)))
        @test length(p.constraints) == length(cl)

        # invalid half-space
        P = Polyhedra.HalfSpace(N[0, 0], N(1)) ∩ Polyhedra.HalfSpace(N[0, 0], N(-1))
        @test_throws AssertionError convert(HPolytope, P)
    end

    # vertices_list of "universal polytope" (strictly speaking: illegal input)
    @test vertices_list(HPolytope{N}()) == Vector{Vector{N}}()

    # 1D vertices_list
    P = convert(HPolytope, Interval(1, 2))
    @test vertices_list(P) == [N[1], N[2]]

    # 2D vertices_list
    P = convert(HPolytope, Hyperrectangle(N[1, 1], N[2, 2]))
    @test vertices_list(P) == [N[3, 3], N[-1, 3], N[-1, -1], N[3, -1]]

    @static if isdefined(@__MODULE__, :Polyhedra)
        # nD (n > 2) vertices_list
        P = convert(HPolytope, Hyperrectangle(N[1, 1, 1], N[2, 2, 2]))
        vlist = vertices_list(P)
        @test ispermutation(vlist,
                            [N[3, 3, -1], N[3, -1, -1], N[-1, -1, -1], N[-1, 3, -1],
                             N[3, 3, 3], N[3, -1, 3], N[-1, -1, 3], N[-1, 3, 3]])
    end

    # isempty
    @test !isempty(HPolytope{N}())  # note: this object is illegal

    # translation
    P = HPolytope([HalfSpace(N[1, 0], N(1)),
                   HalfSpace(N[0, 1], N(1)),
                   HalfSpace(N[-1, -0], N(1)),
                   HalfSpace(N[-0, -1], N(1))])
    P2 = translate(P, N[1, 2])
    @test P2 isa HPolytope{N} && ispermutation(constraints_list(P2),
                                               [HalfSpace(N[1, 0], N(2)), HalfSpace(N[0, 1], N(3)),
                                                HalfSpace(N[-1, -0], N(0)), HalfSpace(N[-0, -1], N(-1))])

    # subset
    H = BallInf(N[0, 0], N(1))
    P = convert(HPolytope, H)
    @test BallInf(N[0, 0], N(1)) ⊆ P
    @test !(BallInf(N[0, 0], N(1.01)) ⊆ P)

    # =====================
    # Concrete linear map
    # =====================
    LM = linear_map(N[2 3; 1 2], P) # invertible matrix
    @test LM isa HPolytope{N}
    @static if isdefined(@__MODULE__, :Polyhedra)
        LM = linear_map(N[2 3; 0 0], P; algorithm="vrep")  # non-invertible matrix
        @test LM isa VPolygon{N}
    end

    M = N[2 1; 0 1]
    L1 = linear_map(M, P; algorithm="inverse")  # calculates inv(M) explicitly
    L2 = linear_map(M, P; algorithm="inverse_right")  # uses transpose(M) \ c.a for each constraint c of P
    L3 = linear_map(M, P; algorithm="vrep")  # uses V-representation
    L3ch = linear_map(M, P; algorithm="vrep_chull")  # uses V-representation and eliminates redundant vertices
    @test_throws ArgumentError linear_map(M, P, algorithm="xyz")  # unknown algorithm
    L4 = linear_map(M, P; cond_tol=1e3)  # set a custom tolerance for the condition number (invertibility check)
    L5 = linear_map(M, P; check_invertibility=false)  # invertibility known
    L6 = linear_map(M, P; inverse=inv(M))  # pass inverse, uses "inverse"
    L7 = linear_map_inverse(inv(M), P)  # convenience function to pass inverse but not M
    p = center(H)
    @test p ∈ P
    @test all([M * p ∈ Li for Li in [L1, L2, L3, L3ch, L4, L5, L6, L7]])
    @test L3 isa VPolygon{N}

    # linear map for mixed types
    M = [2 1; 0 1] # Int's
    LM = linear_map(M, P)
    @test LM isa HPolytope{N}

    if N != Rational{Int}
        M = N[1 1; 2 2; 3 4] # non-invertible matrix of rank 2
        P = HPolytope([HalfSpace(N[1, 0], N(1)),
                       HalfSpace(N[0, 1], N(1)),
                       HalfSpace(N[-1, -0], N(1)),
                       HalfSpace(N[-0, -1], N(1))])
        L7 = linear_map(M, P; algorithm="lift")
        @test L7 isa HPolytope{N}
        @static if isdefined(@__MODULE__, :Polyhedra)
            L7_vrep = linear_map(M, P; algorithm="vrep")
            if N == Float64
                @test L7 ⊆ L7_vrep && L7_vrep ⊆ L7
            end
            # For Float32 we need to support mixed types,
            # ρ(::Vector{Float64}, ::HPolytope{Float32})
        end
    end

    # linear_map_inverse / affine_map_inverse
    M = N[1 0; 0 0]
    # redundant constraint
    P = HPolytope(HalfSpace[HalfSpace(N[1, 1], N(2)), HalfSpace(N[-1, 0], N(0)),
                            HalfSpace(N[0, -1], N(0))])
    X = linear_map_inverse(M, P)
    @test isequivalent(X, HPolyhedron([HalfSpace(N[1, 0], N(2)), HalfSpace(N[-1, 0], N(0))]))
    Y = affine_map_inverse(M, P, zeros(N, 2))
    @test X == Y
    Y = affine_map_inverse(M, P, ones(N, 2))
    @test isequivalent(Y, HPolyhedron([HalfSpace(N[1, 0], N(0)), HalfSpace(N[-1, 0], N(1))]))
    # infeasible constraint
    P = HPolytope(HalfSpace[HalfSpace(N[-1, -1], N(3)), HalfSpace(N[1, 0], N(-1)),
                            HalfSpace(N[0, 1], N(-1))])
    X = linear_map_inverse(M, P)
    @test isempty(X)
    Y = affine_map_inverse(M, P, zeros(N, 2))
    @test X == Y

    # convert
    P = HPolytope([HalfSpace(SingleEntryVector(1, 1, N(1)), N(1)),
                   HalfSpace(SingleEntryVector(1, 1, N(-1)), N(0))])
    Q = convert(typeof(P), P)
    @test P == Q
    @static if isdefined(@__MODULE__, :Polyhedra)
        Q = convert(HPolytope{N,Vector{N}}, tovrep(P))
        @test P == Q
    end
    H = HalfSpace(N[1], N(1))
    @test_throws ErrorException convert(HPolytope, H)
    @test_throws ErrorException convert(HPolytope{N,Vector{N}}, H)
    Z = SimpleSparsePolynomialZonotope(N[2, 0], N[1 2; 2 2.0], [1 4; 1 2])
    @test_throws ErrorException convert(HPolytope, Z)
    @test_throws ErrorException convert(HPolytope{N,Vector{N}}, Z)
    P = convert(HPolytope, EmptySet{N}(2))
    @test P isa HPolytope{N} && dim(P) == 2 && isempty(P)

    # -----
    # V-rep
    # -----

    # constructor from a VPolygon
    polygon = VPolygon([N[0, 0], N[1, 0], N[0, 1]])
    p = convert(VPolytope, polygon)
    @test vertices_list(polygon) == vertices_list(p)

    # dim
    @test dim(p) == 2

    # support function/vector
    d = N[1, 0]
    @test ρ(d, p) == N(1)
    @test σ(d, p) == N[1, 0]
    # empty polytope
    V = VPolytope{N}()
    @test_throws ErrorException ρ(d, V)
    @test_throws ErrorException σ(d, V)
    # one vertex
    V = VPolytope([N[2, 1]])
    @test σ(d, V) == N[2, 1]

    # boundedness
    @test isbounded(p)

    # ispolyhedral
    @test ispolyhedral(p)

    # isempty
    @test !isempty(p)
    @test isempty(VPolytope{N}())

    # vertices_list function
    @test vertices_list(p) == p.vertices

    # low/high/extrema
    p2 = VPolytope([N[5, 1], N[4, 0], N[3, 1], N[4, 2]])
    @test low(p2) == N[3, 0] && low(p2, 1) == N(3) && low(p2, 2) == N(0)
    @test high(p2) == N[5, 2] && high(p2, 1) == N(5) && high(p2, 2) == N(2)
    @test extrema(p2) == (N[3, 0], N[5, 2]) && extrema(p2, 1) == (N(3), N(5)) &&
          extrema(p2, 2) == (N(0), N(2))
    # singleton
    p2 = VPolytope([N[1, 2]])
    @test low(p2) == N[1, 2] && low(p2, 1) == N(1) && low(p2, 2) == N(2)
    @test high(p2) == N[1, 2] && high(p2, 1) == N(1) && high(p2, 2) == N(2)
    @test extrema(p2) == (N[1, 2], N[1, 2]) && extrema(p2, 1) == (N(1), N(1)) &&
          extrema(p2, 2) == (N(2), N(2))
    # empty polytope cannot determine dimension and returns empty vectors
    p2 = VPolytope{N}()
    @test low(p2) == N[]
    @test_throws AssertionError low(p2, 1)
    @test high(p2) == N[]
    @test_throws AssertionError high(p2, 1)
    @test extrema(p2) == (N[], N[])
    @test_throws AssertionError extrema(p2, 1)

    @static if isdefined(@__MODULE__, :Polyhedra)
        # conversion to and from Polyhedra's VRep data structure
        V = convert(VPolytope, polyhedron(p))
        vl = vertices_list(V)
        @test length(p.vertices) == length(vl) && vl ⊆ p.vertices

        # list of constraints of a VPolytope; calculates tohrep
        @test ispermutation(constraints_list(V), constraints_list(tohrep(V)))

        # convert empty VPolytope to a polyhedron
        Vempty = VPolytope()
        @test_throws ErrorException polyhedron(Vempty) # needs to pass the (relative) dim
        Pe = polyhedron(Vempty; relative_dimension=2)

        # remove_redundant_vertices
        V2 = remove_redundant_vertices(V)
        @test V2 == V
    end

    # volume
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test volume(p) == N(1 // 2)
    end

    # translation
    @test translate(p, N[1, 2]) == VPolytope([N[1, 2], N[2, 2], N[1, 3]])
    pp = copy(p)
    @test translate!(pp, N[1, 2]) == VPolytope([N[1, 2], N[2, 2], N[1, 3]]) == pp
    # empty polytope
    V = VPolytope{N}()
    @test translate!(copy(V), N[1]) == V

    # copy (see #1002)
    p, q = [N(1)], [N(2)]
    P = VPolytope([p, q])
    Pcopy = deepcopy(P)
    p[1] = N(5)
    # test that Pcopy is independent of P ( = deepcopy)
    @test Pcopy.vertices[1] == [N(1)]

    # concrete projection
    V = VPolytope([N[0, 0, 1], N[0, 1, 0], N[0, -1, 0], N[1, 0, 0]])
    @test project(V, [1]) == Interval(N(0), N(1))
    @test project(V, [1, 2]) == VPolygon([N[0, 0], N[0, 1], N[0, -1], N[1, 0]])
    @test project(V, [1, 2, 3]) == V
    # empty polytope
    V = VPolytope{N}()
    @test_throws DimensionMismatch project(V, [1])

    # linear_map with redundant vertices
    A = N[1 0; 0 0]
    P = VPolytope([N[1, 1], N[-1, 1], N[1, -1], N[-1, -1]])
    Q1 = linear_map(A, P)
    vlist1 = convex_hull(vertices_list(Q1))
    Q2 = linear_map(A, P; apply_convex_hull=true)
    vlist2 = vertices_list(Q2)
    @test ispermutation(vlist1, [N[1, 0], N[-1, 0]])
    @test ispermutation(vlist1, vlist2)

    # permutation
    P = VPolytope([N[1, 2, 3], N[4, 5, 6]])
    Q = VPolytope([N[3, 1, 2], N[6, 4, 5]])
    @test permute(P, [3, 1, 2]) == Q

    @static if isdefined(@__MODULE__, :StaticArrays)
        # construction with static arrays
        vlist = [SVector{2}(N[0, 0]), SVector{2}(N[1, 0]), SVector{2}(N[0, 1])]
        V = VPolytope(vlist)
        @test ispermutation(vertices_list(V), vlist)
        # check that the outer container can also be static
        vlist = SVector{3}(SVector{2}(N[0, 0]), SVector{2}(N[1, 0]), SVector{2}(N[0, 1]))
        V = VPolytope(vlist)
        @test ispermutation(vertices_list(V), vlist)
    end

    # reflect
    P = VPolytope([N[1, 2, 3], N[4, 5, 6]])
    Q = reflect(P)
    vlist = [N[-1, -2, -3], N[-4, -5, -6]]
    @test Q isa VPolytope{N} && ispermutation(vertices_list(Q), vlist)

    # scale/scale!
    P2 = copy(P)
    scale!(N(2), P2)
    @test scale(N(2), P) == P2 == VPolytope([N[2, 4, 6], N[8, 10, 12]])

    # permute
    P = HPolytope(HalfSpace[HalfSpace(N[1, 2], N(3)), HalfSpace(N[-1, -2], N(-3))])
    @test permute(P, 1:2) == P
    @test permute(P, [2, 1]) ==
          HPolytope(HalfSpace[HalfSpace(N[2, 1], N(3)), HalfSpace(N[-2, -1], N(-3))])
end

# default Float64 constructors
@test HPolytope() isa HPolytope{Float64}
@test VPolytope() isa VPolytope{Float64,Vector{Float64}}

for N in [Float64, Float32]
    # rand
    p = rand(HPolytope; N=N, dim=2, num_vertices=0)
    @test p isa HPolytope{N} && dim(p) == 2 && isempty(p)
    p = rand(HPolytope; N=N, dim=3, num_vertices=0)
    @test p isa HPolytope{N} && dim(p) == 3 && isempty(p)
    p = rand(HPolytope; N=N, num_vertices=1)
    @test p isa HPolytope{N} && dim(p) == 2 && !isempty(p)
    p = rand(HPolytope; N=N, dim=1)
    @test p isa HPolytope{N} && dim(p) == 1
    p = rand(HPolytope; N=N, dim=2)
    @test p isa HPolytope{N} && dim(p) == 2
    @static if isdefined(@__MODULE__, :Polyhedra)
        p = rand(HPolytope; N=N, dim=3)
        @test p isa HPolytope{Float64} && dim(p) == 3
    else
        @test_throws AssertionError rand(HPolytope; N=N, dim=3)
    end
    p = rand(VPolytope; N=N)
    @test p isa VPolytope{N} && dim(p) == 2
    p = rand(VPolytope; N=N, num_vertices=5)
    @test p isa VPolytope{N} && length(p.vertices) <= 5

    # normalization
    p1 = HPolytope([HalfSpace(N[1e5], N(3e5)), HalfSpace(N[-2e5], N(4e5))])
    p2 = normalize(p1)
    for hs in constraints_list(p2)
        @test norm(hs.a) == N(1)
    end
end

for N in [Float64]
    p = HPolytope{N}()
    c1 = LinearConstraint(N[2, 2], N(12))
    c2 = LinearConstraint(N[-3, 3], N(6))
    c3 = LinearConstraint(N[-1, -1], N(0))
    c4 = LinearConstraint(N[2, -4], N(0))
    addconstraint!(p, c3)
    addconstraint!(p, c1)
    addconstraint!(p, c4)
    addconstraint!(p, c2)

    # support vector
    d = N[1, 0]
    @test σ(d, p) == N[4, 2]
    d = N[0, 1]
    @test σ(d, p) == N[2, 4]
    d = N[-1, 0]
    @test σ(d, p) == N[-1, 1]
    d = N[0, -1]
    @test σ(d, p) == N[0, 0]

    # boundedness
    @test isbounded(p) && isbounded(p, false)
    p2 = HPolytope{N}()
    @test isbounded(p2) && !isbounded(p2, false)

    # remove redundant constraints
    P = HPolytope([HalfSpace(N[1, 0], N(1)),
                   HalfSpace(N[0, 1], N(1)),
                   HalfSpace(N[-1, -0], N(1)),
                   HalfSpace(N[-0, -1], N(1)),
                   HalfSpace(N[1, 0], N(2))]) # redundant

    Pred = remove_redundant_constraints(P)
    @test length(Pred.constraints) == 4
    @test length(P.constraints) == 5

    # test in-place removal of redundancies
    remove_redundant_constraints!(P)
    @test length(P.constraints) == 4

    # membership
    p = VPolytope([N[0, 0], N[1, 0], N[0, 1]])
    @test N[0.49, 0.49] ∈ p
    @test N[0.51, 0.51] ∉ p
    q = VPolytope([N[0, 1], N[0, 2]])
    @test N[0, 1 // 2] ∉ q
    # empty polytope
    V = VPolytope{N}()
    @test_throws DimensionMismatch N[0] ∈ V
    # one vertex
    V = VPolytope([N[1, 2]])
    @test N[1, 2] ∈ V && N[2, 2] ∉ V

    # creation from a matrix (each column is a vertex)
    pmat = VPolytope(copy(N[0 0; 1 0; 0 1]'))
    @test p == pmat

    # inclusion (see #1809)
    X = BallInf(N[0.1, 0.2, 0.1], N(0.3))
    Y = convert(HPolytope, X)
    @test Y ⊆ X

    # double inclusion check
    @static if isdefined(@__MODULE__, :Polyhedra)
        Z = convert(VPolytope, Y)
        @test isequivalent(Y, Z)
    end

    # negative double inclusion check
    X_eps = BallInf(N[0.1, 0.2, 0.1], N(0.30001))
    @test !isequivalent(X, X_eps)

    # rectangular map
    M2 = N[2 1; 0 1; 3 3]
    Q = linear_map(M2, P)
    P2 = linear_map_inverse(M2, Q)
    @test isequivalent(P2, P)

    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :CDDLib)
        # -----
        # H-rep
        # -----
        # support function/vector
        d = N[1, 0]
        p_unbounded = HPolytope([LinearConstraint(N[-1, 0], N(0))])
        @test_throws ErrorException σ(d, p_unbounded)
        @test_throws ErrorException ρ(d, p_unbounded)
        p_infeasible = HPolytope([LinearConstraint(N[1], N(0)),
                                  LinearConstraint(N[-1], N(-1))])
        @test_throws ArgumentError σ(N[1], p_infeasible)
        @test_throws ArgumentError ρ(N[1], p_infeasible)

        # empty intersection
        A = [N(0) N(-1); N(-1) N(0); N(1) N(1)]
        b = N[0, 0, 1]
        p1 = HPolytope(A, b)
        A = [N(0) N(1); N(1) N(0); N(-1) N(-1)]
        b = N[2, 2, -2]
        p2 = HPolytope(A, b)
        @test intersection(p1, p2) isa EmptySet{N}
        @test intersection(p1, p2; backend=LazySets.default_polyhedra_backend(p1)) isa EmptySet{N}
        @test intersection(p1, p2; backend=CDDLib.Library()) isa EmptySet{N}

        # intersection with half-space
        hs = HalfSpace(N[2, -2], N(-1))
        c1 = constraints_list(intersection(p1, hs))
        c2 = constraints_list(intersection(hs, p1))
        @test length(c1) == 3 && ispermutation(c1, c2)

        # Cartesian product
        A = [N(1) N(-1)]'
        b = N[1, 0]
        p1 = HPolytope(A, b)
        p2 = copy(p1)
        cp = cartesian_product(p1, p2)
        cl = constraints_list(cp)
        @test length(cl) == 4

        # vertices_list
        vl = vertices_list(p1)
        @test ispermutation(vl, [N[0], N[1]])

        # isempty
        @test !isempty(p)
        P = HPolytope([LinearConstraint(N[1, 0], N(0))])      # x <= 0
        @test !isempty(P)
        addconstraint!(P, LinearConstraint(N[-1, 0], N(-1)))  # x >= 1
        @test isempty(P)

        # checking for empty intersection (also test symmetric calls)
        P = convert(HPolytope, BallInf(zeros(N, 2), N(1)))
        Q = convert(HPolytope, BallInf(ones(N, 2), N(1)))
        R = convert(HPolytope, BallInf(3 * ones(N, 2), N(1)))
        res, w = isdisjoint(P, Q, true)
        @test !isdisjoint(P, Q) && !res && w ∈ P && w ∈ Q
        @test !isdisjoint(P, Q; algorithm="sufficient")
        res, w = isdisjoint(Q, P, true)
        @test !isdisjoint(Q, P) && !res && w ∈ P && w ∈ Q
        res, w = isdisjoint(Q, R, true)
        @test !isdisjoint(Q, R) && !res && w ∈ Q && w ∈ R
        res, w = isdisjoint(R, Q, true)
        @test !isdisjoint(R, Q) && !res && w ∈ Q && w ∈ R
        res, w = isdisjoint(P, R, true)
        @test isdisjoint(P, R) && res && w == N[]
        res, w = isdisjoint(R, P, true)
        @test isdisjoint(R, P) && res && w == N[]
        res, w = isdisjoint(P, R, true; algorithm="sufficient")
        @test isdisjoint(P, R; algorithm="sufficient") && res && w == N[]

        # test that one can pass a sparse vector as the direction (see #1011)
        P = HPolytope([HalfSpace(N[1, 0], N(1)),
                       HalfSpace(N[0, 1], N(1)),
                       HalfSpace(N[-1, -1], N(-1))])
        @test an_element(P) ∈ P

        # tovrep from HPolytope
        A = N[0 1; 1 0; -1 -1]
        b = N[0.25, 0.25, -0]
        P = HPolytope(A, b)
        @test tovrep(P) isa VPolytope{N}
        @test tohrep(P) isa HPolytope{N}  # no-op

        # linear map
        H = BallInf(N[0, 0], N(1))
        P = convert(HPolytope, H)
        M = N[2 3; 1 2]
        MP = linear_map(M, P)  # check concrete linear map
        @test M * an_element(P) ∈ MP   # check that the image of an element is in the linear map

        # check boundedness after conversion
        H = Hyperrectangle(N[1, 1], N[2, 2])
        HPolytope(constraints_list(H); check_boundedness=true)

        # -----
        # V-rep
        # -----

        # intersection
        p1 = VPolytope(vertices_list(BallInf(N[0, 0], N(1))))
        p2 = VPolytope(vertices_list(BallInf(N[1, 1], N(1))))
        cap = intersection(p1, p2)
        vlist = [N[1, 1], N[0, 1], N[1, 0], N[0, 0]]
        @test ispermutation(vertices_list(cap), vlist)
        # other polytopic sets
        p3 = VPolygon(vertices_list(p2))
        cap = intersection(p1, p3)
        @test ispermutation(vertices_list(cap), vlist)

        # 2D intersection
        paux = VPolytope([N[0, 0], N[1, 0], N[0, 1], N[1, 1]])
        qaux = VPolytope([N[1, -1 / 2], N[-1 / 2, 1], N[-1 / 2, -1 / 2]])
        xaux = intersection(paux, qaux)
        oaux = VPolytope([N[0, 0], N[1 / 2, 0], N[0, 1 / 2]])
        @test xaux ⊆ oaux && oaux ⊆ xaux # TODO use isequivalent

        # mixed types
        paux = VPolygon([N[0, 0], N[1, 0], N[0, 1], N[1, 1]])
        qaux = VPolytope([N[1, -1 / 2], N[-1 / 2, 1], N[-1 / 2, -1 / 2]])
        xaux = intersection(paux, qaux)
        oaux = VPolytope([N[0, 0], N[1 / 2, 0], N[0, 1 / 2]])
        @test xaux ⊆ oaux && oaux ⊆ xaux # TODO use isequivalent

        # 1D set
        paux = VPolytope([N[0], N[1]])
        qaux = VPolytope([N[-1 / 2], N[1 / 2]])
        xaux = intersection(paux, qaux)
        oaux = VPolytope([N[0], N[1 / 2]])
        @test xaux ⊆ oaux && oaux ⊆ xaux # TODO use isequivalent

        # isuniversal
        answer, w = isuniversal(p1, true)
        @test !isuniversal(p1) && !answer && w ∉ p1

        # convex hull
        v1 = N[1, 0]
        v2 = N[0, 1]
        v3 = N[-1, 0]
        v4 = N[0, -1]
        v5 = N[0, 0]
        p1 = VPolytope([v1, v2, v5])
        p2 = VPolytope([v3, v4])
        ch = convex_hull(p1, p2)
        vl = vertices_list(ch)
        @test ispermutation(vl, [v1, v2, v3, v4])

        # Cartesian product
        p1 = VPolytope([N[0, 0], N[1, 1]])
        p2 = VPolytope([N[2]])
        cp = cartesian_product(p1, p2)
        vl = vertices_list(cp)
        @test ispermutation(vl, [N[0, 0, 2], N[1, 1, 2]])

        # tohrep from VPolytope
        P = VPolytope([v1, v2, v3, v4, v5])
        @test tohrep(P) isa HPolytope{N}
        @test tovrep(P) isa VPolytope{N}  # no-op

        # subset (see #974)
        H = BallInf(N[0, 0], N(1))
        P = convert(VPolytope, H)
        @test BallInf(N[0, 0], N(1)) ⊆ P
        @test !(BallInf(N[0, 0], N(1.01)) ⊆ P)

        # concrete minkowski sum
        for T in (HPolytope, VPolytope)
            B = convert(T, BallInf(N[0, 0, 0], N(1)))
            X = minkowski_sum(B, B)
            twoB = 2.0 * B
            @test X ⊆ twoB && twoB ⊆ X
        end

        # concrete Reflection
        F4 = N[3 0; 0 7]
        g4 = N[1, 2]
        P4 = HPolyhedron(F4, g4)
        F, g = tosimplehrep(reflect(reflect(P4)))
        @test F4 == F && g4 == g

        # same but specifying a custom polyhedral computations backend (CDDLib)
        X = minkowski_sum(B, B; backend=CDDLib.Library())
        twoB = 2.0 * B
        @test X ⊆ twoB && twoB ⊆ X

        P1 = VPolytope([N[0, 0, 0], N[0, 1, 0]])
        P2 = VPolytope([N[0, 0, 0], N[1, 0, 0]])
        Q = minkowski_sum(P1, P2)
        @test ispermutation(vertices_list(Q), [N[0, 0, 0], N[0, 1, 0], N[1, 0, 0], N[1, 1, 0]])

        # fallback conversion to vertex representation
        B3 = BallInf(zeros(N, 3), N(1))
        U = Matrix(N(1) * I, 3, 3) * B3
        Uv = convert(VPolytope, U; prune=false)
        @test ispermutation(vertices_list(Uv), vertices_list(B3))

        # -----------------
        # mixed H-rep/V-rep
        # -----------------

        # intersection
        p1 = convert(HPolytope, BallInf(N[0, 0], N(1)))
        p2 = convert(VPolytope, BallInf(N[1, 1], N(1)))
        cap1 = intersection(p1, p2)
        cap2 = intersection(p2, p1)
        @test ispermutation([round.(v) for v in vertices_list(cap1)],
                            [N[1, 1], N[0, 1], N[1, 0], N[0, 0]])
        @test ispermutation([round.(v) for v in vertices_list(cap1)],
                            [round.(v) for v in vertices_list(cap2)])

        # Chebyshev center
        B = BallInf(N[0, 0], N(1))  # Chebyshev center is unique
        c1, r1 = chebyshev_center_radius(B)
        P = convert(HPolytope, B)
        c2, r2 = chebyshev_center_radius(P)
        @test c1 == c2 == center(B) && c1 isa AbstractVector{N}
        @test r1 == r2 == B.radius

        # concrete projection
        πP = project(P, [1])
        @test πP isa HPolytope{N}
        @test ispermutation(constraints_list(πP), [HalfSpace(N[-1], N(1)),
                                                   HalfSpace(N[1], N(1))])
    end

    # tests that require Symbolics
    @static if isdefined(@__MODULE__, :Symbolics)
        vars = @variables x y
        p1 = HPolytope([x + y <= 1, x + y >= -1, x - y <= 1, x - y >= -1], vars)
        p2 = HPolytope([x + y <= 1, x + y >= -1, x - y <= 1, x - y >= -1])
        @test p1 == p2
        p2 = HPolytope([x + y <= 1, x + y >= -1, x - y <= 1, x - y >= -1], [vars])
        @test p1 == p2
        b1 = Ball1(zeros(2), 1.0)
        @test isequivalent(p1, b1)
    end

    # concrete projection of a polytope (see issue #2536)
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :CDDLib)
        if N == Float64
            X = HPolytope([HalfSpace([1.0, 1.0, 0.0], 4.0),
                           HalfSpace([-1.0, -1.0, -0.0], -4.0),
                           HalfSpace([-1.0, 0.0, 0.0], -0.0),
                           HalfSpace([0.0, -1.0, 0.0], -0.0),
                           HalfSpace([0.0, 0.0, -1.0], -0.0),
                           HalfSpace([1.0, 0.0, 0.0], 10.0),
                           HalfSpace([0.0, 1.0, 0.0], 10.0),
                           HalfSpace([0.0, 0.0, 1.0], 10.0),
                           HalfSpace([0.0, 1.0, 1.0], 6.0)])
            v12 = [N[0, 4], N[4, 0]]
            @test ispermutation(vertices_list(project(X, 1:2)), v12)
            @test ispermutation(vertices_list(overapproximate(Projection(X, [1, 2]), 1e-3)), v12)
        end
    end
end

@test !isoperationtype(HPolytope)
@test !isoperationtype(VPolytope)

for N in [Float64, Rational{Int}]
    @static if isdefined(@__MODULE__, :MiniQhull)
        # Delaunay triangulation
        vlist = [N[0, 0, 0], N[0, 0, 1], N[0, 1, 0], N[1, 0, 0]]  # tetrahedron
        V = VPolytope(vlist)
        D = triangulate(V; algorithm="delaunay")
        @test length(D) == 1 && isequivalent(array(D)[1], V)
        D = triangulate(V; algorithm="delaunay", compute_triangles_3d=true)
        @test length(D) == 4
        for P in array(D)
            @test isequivalent(P, VPolytope(vlist[[1, 2, 3]])) ||
                  isequivalent(P, VPolytope(vlist[[1, 2, 4]])) ||
                  isequivalent(P, VPolytope(vlist[[1, 3, 4]])) ||
                  isequivalent(P, VPolytope(vlist[[2, 3, 4]]))
        end
    end
end
