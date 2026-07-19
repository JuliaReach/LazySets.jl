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
    return _cartesian_product_vrep(P1, P2; backend1=backend, backend2=backend)
end

# see ext/LazySetsPolyhedraExt.jl
function _cartesian_product_vrep(P1, P2; backend1=nothing, backend2=nothing)
    mod = Base.get_extension(@__MODULE__, :LazySetsPolyhedraExt)
    require(mod, :Polyhedra; fun_name="cartesian_product")
    error()
end
