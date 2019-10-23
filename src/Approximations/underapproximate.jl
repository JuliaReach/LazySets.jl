"""
    underapproximate(X::LazySet{N}, dirs::AbstractDirections;
                    [apply_convex_hull]::Bool=false) where {N<:Real}

Compute the underapproximation of a convex set by sampling support vectors.

### Input

- `X`                 -- set
- `dirs`              -- directions
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the support vectors with a convex hull operation
### Output

The `VPolytope` obtained by taking the convex hull of the support vectors of `X`
along the directions determined by `dirs`.
"""
function underapproximate(X::LazySet{N}, dirs::AbstractDirections;
                          apply_convex_hull::Bool=false) where {N<:Real}
    @assert LazySets.dim(X) == LazySets.dim(dirs)
    vinner = Vector{Vector{N}}(undef, length(dirs)) # TODO could be generalized to eltype(dirs)
    @inbounds for (i, di) in enumerate(dirs)
        vinner[i] = σ(Vector(di), X)
    end
    if apply_convex_hull
        convex_hull!(vinner)
    end
    return VPolytope(vinner)
end

"""
    underapproximate(P::AbstractPolyhedron{N}, ::Type{Ellipsoid};
                     backend=default_SDP_solver(N),
                     interior_point::AbstractVector{N}=chebyshev_center(P)
                    ) where {N<:Real}

Underapproximate a polyhedron by the maximum-volume ellipsoid.

### Input

- `P`         -- polyhedral set
- `Ellipsoid` -- type for dispatch
- `backend`   -- (optional, default: `default_SDP_solver(N)`) backend to solve
                 the semi-definite program
- `interior_point` -- (optional, default: `nothing`) an interior point of the
                      ellipsoid (needed for the algorithm: see below for
                      details)

### Output

An ellipsoid ``E`` such that ``E ⊆ P`` and ``E`` has maximal volume.

### Notes

The maximum-volume ellipsoid is also called *Löwner-John ellipsoid*.

The algorithm requires to specify a point in the interior of the ellipsoid which
can be specified with the argument `interior_point`.
If no such point is known (i.e., `interior_point == nothing`), the
implementation will compute the Chebyshev center (see
[`chebyshev_center`](@ref)).

The Chebyshev center will be computed using the default LP solver because it is
not guaranteed that the SDP solver `backend` supports LPs.
If a custom LP solver should be used, the Chebyshev center needs to be computed
and passed manually.

### Algorithm

TODO
"""
function underapproximate(P::AbstractPolyhedron{N}, ::Type{Ellipsoid};
                          backend=default_SDP_solver(N),
                          interior_point::AbstractVector{N}=chebyshev_center(P)
                         ) where {N<:Real}
    require(:SetProg; fun_name="underapproximate")
end

function load_setprog_underapproximate()
return quote

function underapproximate(P::AbstractPolyhedron{N}, ::Type{Ellipsoid};
                          backend=default_SDP_solver(N),
                          interior_point::AbstractVector{N}=chebyshev_center(P)
                         ) where {N<:Real}
    # create SDP model
    model = SDPModel{N}()

    # create symbolic ellipsoid S
    SetProg.@variable(model, S, Ellipsoid(point=SetProg.InteriorPoint(interior_point)))

    # add subset constraint
    Q = polyhedron(P)  # convert P to a Polyhedra.polyhedron
    SetProg.@constraint(model, S ⊆ Q)

    # add maximum-volume objective
    SetProg.@objective(model, Max, nth_root(volume(S)))

    # solve SDP
    MOI.copy_to(backend, model)  # send model to solver
    MOI.optimize!(backend)
    ellipsoid = SetProg.value(S)

    # convert result to a LazySets.Ellipsoid
    return convert(Ellipsoid, ellipsoid)
end

end # quote
end # load_setprog_underapproximate
