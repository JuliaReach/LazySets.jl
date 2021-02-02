export cartesian_product

"""
    cartesian_product(X::LazySet, Y::LazySet; [backend]=nothing, [algorithm]::String="vrep")

Compute the Cartesian product of two sets.

### Input

- `X`         -- set
- `Y`         -- another set
- `backend`   -- (optional, default: `nothing`) the polyhedral computations backend
- `algorithm` -- (optional, default: "hrep") the method used to transform each set
                 `X` and `Y` before taking the Cartesian product; choose between
                 "vrep" (use the vertex representation) and "hrep" (use the constraint representation)

### Output

The `VPolytope` (if "vrep" was used) or `HPolytope` (if "hrep" was used)
obtained by the concrete Cartesian product of `X` and `Y`.

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
        Xp = HPolyhedron(constraints_list(X))
        Yp = HPolyhedron(constraints_list(Y))
        Pout = cartesian_product(Xp, Yp, backend=backend)

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

"""
    cartesian_product(P1::HPoly, P2::HPoly; [backend]=nothing)

Compute the Cartesian product of two polyhedra in H-representaion.

### Input

- `P1`         -- polyhedron
- `P2`         -- another polyhedron
- `backend`    -- (optional, default: `nothing`) the polyhedral computations backend

### Output

The polyhedron obtained by the concrete cartesian product of `P1` and `P2`.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/).
"""
function cartesian_product(P1::HPoly, P2::HPoly; backend=nothing)
    require(:Polyhedra; fun_name="`cartesian_product")

    return _cartesian_product_hrep(P1, P2, backend1=backend, backend2=backend)
end

function _cartesian_product_hrep(P1, P2; backend1, backend2)
    if isnothing(backend1)
        backend1 = default_polyhedra_backend(P1)
    end
    if isnothing(backend2)
        backend2 = default_polyhedra_backend(P2)
    end

    P1′ = polyhedron(P1; backend=backend1)
    P2′ = polyhedron(P2; backend=backend2)
    Pout = Polyhedra.hcartesianproduct(P1′, P2′)
    return convert(basetype(P1), Pout)
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
