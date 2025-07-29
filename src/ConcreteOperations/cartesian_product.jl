"""
    cartesian_product(X::LazySet, Y::LazySet; [backend]=nothing,
                      [algorithm]::String="vrep")

Compute the Cartesian product of two sets.

### Input

- `X`         -- set
- `Y`         -- set
- `backend`   -- (optional, default: `nothing`) backend for polyhedral
                 computation
- `algorithm` -- (optional, default: "hrep") the method used to transform the
                 sets `X` and `Y` before taking the Cartesian product; choose
                 between:

    - "vrep" (use the vertex representation)
    - "hrep" (use the constraint representation)
    - "hrep_polyhedra" (use the constraint representation and `Polyhedra`)

### Output

The `VPolytope` (if `"vrep"` was used) or `HPolytope`/`HPolyhedron` (if `"hrep"`
or `"hrep_polyhedra"` was used) obtained by the concrete Cartesian product of
`X` and `Y`. The choice between `HPolytope` and `HPolyhedron` is made based on
boundedness information deduced from the type.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/).

If `X` can be converted to a one-dimensional interval and the vertices of `Y`
are available, use `algorithm="vrep"`.
"""
function cartesian_product(X::LazySet, Y::LazySet; backend=nothing,
                           algorithm::String="hrep")
    if algorithm == "vrep"
        Yv = VPolytope(vertices_list(Y))
        if dim(X) == 1 && isconvextype(typeof(X))  # special case
            Xv = convert(Interval, X)  # works because X is convex
            return cartesian_product(Xv, Yv)
        end
        Xv = VPolytope(vertices_list(X))
        Pout = cartesian_product(Xv, Yv; backend=backend)

    elseif algorithm == "hrep"
        Pout = _cartesian_product_hrep(X, Y)

    elseif algorithm == "hrep_polyhedra"
        Xp = HPolyhedron(constraints_list(X))
        if isempty(Xp.constraints) && isuniversal(X)  # `isuniversal` needed for empty HPolytope
            Xp = X isa Universe ? X : Universe(dim(X))
        end
        Yp = HPolyhedron(constraints_list(Y))
        if isempty(Yp.constraints) && isuniversal(Y)  # `isuniversal` needed for empty HPolytope
            Yp = Y isa Universe ? Y : Universe(dim(Y))
        end
        Pout = _cartesian_product_hrep_polyhedra(Xp, Yp; backend1=backend,
                                                 backend2=backend)

    else
        throw(ArgumentError("expected algorithm `\"vrep\"`, `\"hrep\"`, or " *
                            "`\"hrep_polyhedra\"`, but got `$algorithm`"))
    end
    return Pout
end

function _cartesian_product_hrep(X::S1, Y::S2) where {S1<:LazySet,S2<:LazySet}
    N = promote_type(eltype(X), eltype(Y))
    U1 = Universe{N}(dim(X))
    clist1 = [cartesian_product(U1, c) for c in constraints_list(Y)]
    U2 = Universe{N}(dim(Y))
    clist2 = [cartesian_product(c, U2) for c in constraints_list(X)]
    clist = vcat(clist1, clist2)
    if isboundedtype(S1) && isboundedtype(S2)
        return HPolytope(clist)
    else
        return HPolyhedron(clist)
    end
end

function _cartesian_product_hrep_polyhedra(P1::PT1, P2::PT2; backend1=nothing,
                                           backend2=nothing) where {PT1,PT2}
    require(@__MODULE__, :Polyhedra; fun_name="`cartesian_product")

    if isnothing(backend1)
        backend1 = default_polyhedra_backend(P1)
    end
    if isnothing(backend2)
        backend2 = default_polyhedra_backend(P2)
    end

    P1′ = polyhedron(P1; backend=backend1)
    P2′ = polyhedron(P2; backend=backend2)
    Pout = Polyhedra.hcartesianproduct(P1′, P2′)

    PT = isboundedtype(PT1) && isboundedtype(PT2) ? HPolytope : HPolyhedron
    return convert(PT, Pout)
end

# if the first set is an interval => the result is always a polytope
function cartesian_product(I::Interval, P::Union{VPolygon,VPolytope})
    a = min(I)
    b = max(I)
    vlist = vertices_list(P)
    m = length(vlist)

    vout = Vector{eltype(vlist)}(undef, 2m)
    @inbounds for i in 1:m
        vout[i] = vcat(a, vlist[i])
        vout[i + m] = vcat(b, vlist[i])
    end
    return VPolytope(vout)
end

function cartesian_product(S1::AbstractSingleton, S2::AbstractSingleton)
    v1 = element(S1)
    v2 = element(S2)
    return Singleton(vcat(v1, v2))
end

function cartesian_product(Z1::AbstractZonotope, Z2::AbstractZonotope)
    N = promote_type(eltype(Z1), eltype(Z2))
    z1 = zeros(N, dim(Z1), ngens(Z2))
    z2 = zeros(N, dim(Z2), ngens(Z1))
    g1 = genmat(Z1)
    g2 = genmat(Z2)
    c = vcat(center(Z1), center(Z2))
    G = [g1 z1; z2 g2]
    return Zonotope(c, G)
end

function cartesian_product(H1::AbstractHyperrectangle,
                           H2::AbstractHyperrectangle)
    c = vcat(center(H1), center(H2))
    r = vcat(radius_hyperrectangle(H1), radius_hyperrectangle(H2))
    return Hyperrectangle(c, r)
end

function cartesian_product(P::AbstractPolyhedron, U::Universe)
    clist = [cartesian_product(H, U) for H in constraints_list(P)]
    return HPolyhedron(clist)
end

function cartesian_product(U::Universe, P::AbstractPolyhedron)
    clist = [cartesian_product(U, H) for H in constraints_list(P)]
    return HPolyhedron(clist)
end

function cartesian_product(H::HalfSpace, U::Universe)
    return HalfSpace(append_zeros(H.a, dim(U)), H.b)
end

function cartesian_product(U::Universe, H::HalfSpace)
    return HalfSpace(prepend_zeros(H.a, dim(U)), H.b)
end

@commutative function cartesian_product(∅::EmptySet, X::LazySet)
    return _cartesian_product_emptyset(∅, X)
end

"""
# Extended help

    cartesian_product(SPZ::AbstractSparsePolynomialZonotope, Z::AbstractZonotope)

### Algorithm

This method implements [Kochdumper21a; Proposition 3.1.22](@citet).
"""
@commutative function cartesian_product(SPZ::AbstractSparsePolynomialZonotope, Z::AbstractZonotope)
    c = vcat(center(SPZ), center(Z))
    G1 = genmat_dep(SPZ)
    G = vcat(G1, zeros(eltype(G1), size(G1)))
    GI = cat(genmat_indep(SPZ), genmat(Z); dims=(1, 2))
    E = expmat(SPZ)
    return SparsePolynomialZonotope(c, G, GI, E, _indexvector(SPZ))
end

"""
# Extended help

    cartesian_product(P1::AbstractSparsePolynomialZonotope, P2::AbstractSparsePolynomialZonotope)

### Algorithm

This method implements [Kochdumper21a; Proposition 3.1.22](@citet).
"""
function cartesian_product(P1::AbstractSparsePolynomialZonotope,
                           P2::AbstractSparsePolynomialZonotope)
    c = vcat(center(P1), center(P2))
    G = cat(genmat_dep(P1), genmat_dep(P2); dims=(1, 2))
    GI = cat(genmat_indep(P1), genmat_indep(P2); dims=(1, 2))
    E = cat(expmat(P1), expmat(P2); dims=(1, 2))
    return SparsePolynomialZonotope(c, G, GI, E)
end
