using Test, LazySets
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # approximate rectification by rectifying the vertices
    # - exact approximation
    X = Ball1(N[1, 1], N(1))
    Y = approximate(Rectification(X))
    @test isequivalent(Y, X)
    # - underapproximation
    X = Ball1(N[2, 2], N(3))
    Y = approximate(Rectification(X))
    Z = VPolygon([N[5, 2], N[2, 5], N[2, 0], N[0, 2]])
    @test Y ⊆ Z && Z ⊆ Y
    # - overapproximation
    X = Ball1(N[-1, -1], N(2))
    Y = approximate(Rectification(X))
    Z = VPolygon([N[0, 0], N[0, 1], N[1, 0]])
    @test Y ⊆ Z && Z ⊆ Y
    # - neither over- nor underapproximation
    X = VPolygon([N[3, 2], N[1, -1], N[-1, 2], N[-1, 4]])
    Y = approximate(Rectification(X))
    Z = VPolygon([N[3, 2], N[1, 0], N[0, 2], N[0, 4]])
    @test Y ⊆ Z && Z ⊆ Y
end
