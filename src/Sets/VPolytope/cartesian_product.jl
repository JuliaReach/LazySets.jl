"""
# Extended help

    cartesian_product(P1::VPolytope, P2::VPolytope; [backend]=nothing)

### Input

- `backend` -- (optional, default: `nothing`) backend for polyhedral computation

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/).
"""
function cartesian_product(P1::VPolytope, P2::VPolytope; backend=nothing)
    require(@__MODULE__, :Polyhedra; fun_name="cartesian_product")

    return _cartesian_product_vrep(P1, P2; backend1=backend, backend2=backend)
end

function _cartesian_product_vrep(P1, P2; backend1=nothing, backend2=nothing)
    if isnothing(backend1)
        backend1 = default_polyhedra_backend(P1)
    end
    if isnothing(backend2)
        backend2 = default_polyhedra_backend(P2)
    end

    P1′ = LazySets.polyhedron(P1; backend=backend1)
    P2′ = LazySets.polyhedron(P2; backend=backend2)
    Pout = Polyhedra.vcartesianproduct(P1′, P2′)
    return convert(VPolytope, Pout)
end
