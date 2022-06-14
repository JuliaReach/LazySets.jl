export cartesian_product

"""
    cartesian_product(X::LazySet, Y::LazySet; [backend]=nothing, [algorithm]::String="vrep")

Compute the Cartesian product of two sets.

### Input

- `X`         -- set
- `Y`         -- another set
- `backend`   -- (optional, default: `nothing`) the polyhedral computations backend
- `algorithm` -- (optional, default: "hrep") the method used to transform each set
                 `X` and `Y` before taking the Cartesian product; choose between:
                 - "vrep" (use the vertex representation)
                 - "hrep" (use the constraint representation)
                 - "hrep_polyhedra" (use the constraint representation and
                   `Polyhedra`)

### Output

The `VPolytope` (if `"vrep"` was used) or `HPolytope` (if `"hrep"` or
`"hrep_polyhedra"` was used) obtained by the concrete Cartesian product of `X`
and `Y`.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/).

If `X` can be converted to a one-dimensional interval and the vertices of `Y`
are available use `algorithm="vrep"`.
"""
function cartesian_product(X::LazySet, Y::LazySet; backend=nothing, algorithm::String="hrep")

    if algorithm == "vrep"
        Yv = VPolytope(vertices_list(Y))
        if dim(X) == 1 # special case
            Xv = convert(Interval, X)  # works if X is convex
            return cartesian_product(Xv, Yv)
        end
        Xv = VPolytope(vertices_list(X))
        Pout = cartesian_product(Xv, Yv, backend=backend)

    elseif algorithm == "hrep"
        Pout = _cartesian_product_hrep(X, Y)

    elseif algorithm == "hrep_polyhedra"
        Xp = HPolyhedron(constraints_list(X))
        if isempty(Xp.constraints) && isuniversal(X)
            Xp = X isa Universe ? X : Universe(dim(X))
        end
        Yp = HPolyhedron(constraints_list(Y))
        if isempty(Yp.constraints) && isuniversal(Y)
            Yp = Y isa Universe ? Y : Universe(dim(Y))
        end
        Pout = _cartesian_product_hrep_polyhedra(Xp, Yp; backend1=backend,
                                                         backend2=backend)

    else
        throw(ArgumentError("expected algorithm `vrep` or `hrep`, got $algorithm"))
    end
    return Pout
end

"""
    cartesian_product(P1::VPolytope, P2::VPolytope; [backend]=nothing)

Compute the Cartesian product of two polytopes in V-representation.

### Input

- `P1`      -- polytope
- `P2`      -- another polytope
- `backend` -- (optional, default: `nothing`) the polyhedral computations backend

### Output

The `VPolytope` obtained by the concrete Cartesian product of `P1` and `P2`.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/).
"""
function cartesian_product(P1::VPolytope, P2::VPolytope; backend=nothing)
    require(:Polyhedra; fun_name="cartesian_product")

    return _cartesian_product_vrep(P1, P2, backend1=backend, backend2=backend)
end

function _cartesian_product_vrep(P1, P2;
                                 backend1=default_polyhedra_backend(P1),
                                 backend2=default_polyhedra_backend(P2))
    if isnothing(backend1)
        backend1 = default_polyhedra_backend(P1)
    end
    if isnothing(backend2)
        backend2 = default_polyhedra_backend(P2)
    end

    P1′ = polyhedron(P1; backend=backend1)
    P2′ = polyhedron(P2; backend=backend2)
    Pout = Polyhedra.vcartesianproduct(P1′, P2′)
    return VPolytope(Pout)
end

function cartesian_product(U1::Universe, U2::Universe)
    return Universe(dim(U1) + dim(U2))
end

function _cartesian_product_hrep(X::S1, Y::S2) where {S1<:LazySet{N}, S2<:LazySet} where {N}
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

function _cartesian_product_hrep_polyhedra(P1::PT1, P2::PT2; backend1, backend2) where {PT1, PT2}
    require(:Polyhedra; fun_name="`cartesian_product")

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
function cartesian_product(I::Interval, P::Union{VPolygon, VPolytope})
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

function cartesian_product(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle)
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

"""
    cartesian_product(P1::SimpleSparsePolynomialZonotope, P2::SimpleSparsePolynomialZonotope)

computes the cartesian product of the simple sparse polynomial zonotopes `P1` and `P2`.

### Input

- `P1` -- simple sparse polynomial zonotope
- `P2` -- simple sparse polynomial zonotope

### Output

The cartesian product of `P1` and `P2`
"""
function cartesian_product(P1::SimpleSparsePolynomialZonotope, P2::SimpleSparsePolynomialZonotope)
    c = vcat(center(P1), center(P2))
    G = cat(genmat(P1), genmat(P2), dims=(1, 2))
    E = cat(expmat(P1), expmat(P2), dims=(1, 2))
    return SimpleSparsePolynomialZonotope(c, G, E)
end
