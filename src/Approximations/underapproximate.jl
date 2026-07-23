"""
    underapproximate(X::LazySet{N}, dirs::AbstractDirections;
                    [apply_convex_hull]::Bool=false) where {N}

Compute the underapproximation of a convex set by sampling support vectors.

### Input

- `X`                 -- convex set
- `dirs`              -- directions
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the support vectors with a convex hull operation
### Output

The `VPolytope` obtained by taking the convex hull of support vectors of `X`
along the directions determined by `dirs`.

### Notes

Since the support vectors are not always unique, this algorithm may return
a strict underapproximation even if the set can be exactly approximated using
the given template.
"""
function underapproximate(X::LazySet{N}, dirs::AbstractDirections;
                          apply_convex_hull::Bool=false) where {N}
    @assert isconvex(X) "this underapproximation requires a convex set"
    @assert dim(X) == dim(dirs) "the dimension of the set ($(dim(X))) does not match the " *
                                "dimension of the directions ($(dim(dirs)))"

    vlist = Vector{Vector{N}}(undef, length(dirs))
    j = 0
    @inbounds for di in dirs
        v = σ(di, X)
        if !any(isinf, v)  # skip directions where the set is unbounded
            j += 1
            vlist[j] = v
        end
    end
    resize!(vlist, j)  # remove entries that were skipped

    if apply_convex_hull
        convex_hull!(vlist)
    end
    return VPolytope(vlist)
end

function underapproximate(X::LazySet, dirs::Type{<:AbstractDirections}; kwargs...)
    return underapproximate(X, dirs(dim(X)), kwargs...)
end

"""
    underapproximate(X::LazySet, ::Type{<:Hyperrectangle};
                     solver=nothing) where {N}

Underapproximate a convex polygon with a hyperrectangle of maximal area.

### Input

- `X`              -- convex polygon
- `Hyperrectangle` -- type for dispatch
- `solver`         -- (optional; default: `nothing`) nonlinear solver; if
                      `nothing`, `default_nln_solver(N)` will be used

### Output

A hyperrectangle underapproximation with maximal area.

### Notes

The implementation only works for 2D sets, but the algorithm can be generalized.

Due to numerical issues, the result may be slightly outside the set.

### Algorithm

The algorithm is taken from [Behroozi19; Theorem 17](@citet) and solves a convex program (in fact
a linear program with nonlinear objective). (The objective is modified to an
equivalent version due to solver issues.)
"""
function underapproximate(X::LazySet{N}, ::Type{<:Hyperrectangle};
                          solver=nothing) where {N}
    if isnothing(solver)
        solver = default_nln_solver(N)
    end
    return _underapproximate_box(X, solver)
end

function _underapproximate_box(X::S, solver) where {S<:LazySet}
    if !ispolytopic(X)
        throw(ArgumentError("box underapproximation is only available for convex polytopes"))
    end

    n = dim(X)
    @assert n == 2 "currently only 2D sets are supported"
    P, b = tosimplehrep(X)

    model = Model(solver)
    set_silent(model)
    @variable(model, u[1:n])
    @variable(model, v[1:n])
    @variable(model, x[1:n])
    @variable(model, y[1:n])
    @variable(model, z[1:n])

    # constraints
    @constraint(model, u[2] == 0)
    @constraint(model, v[1] == 0)
    @constraint(model, [i in 1:n], u[i] - y[i] + x[i] == 0)
    @constraint(model, [i in 1:n], v[i] - z[i] + x[i] == 0)
    @constraint(model, P * x .≤ b)
    @constraint(model, P * y .≤ b)
    @constraint(model, P * z .≤ b)
    @constraint(model, P * (y + z - x) .≤ b)

    # maximize area
    @NLobjective(model, Max, u[1] * v[2])
    optimize!(model)

    # convert solution to hyperrectangle
    @inbounds begin
        l = value.(x)
        h = [value(y[1]), value(z[2])]
        if l[1] > h[1]  # sometimes l > h, then just swap them
            tmp = l
            l = h
            h = tmp
        end
    end
    return Hyperrectangle(; low=l, high=h)
end

"""
    underapproximate(X::LazySet, ::Type{<:Ball2})

Compute the largest inscribed Euclidean ball in a set `X`.

### Input

- `X`     -- set
- `Ball2` -- target type

### Output

A largest `Ball2` contained in `X`.

### Algorithm

We use `chebyshev_center_radius(X)`.
"""
function underapproximate(X::LazySet, ::Type{<:Ball2})
    c, r = chebyshev_center_radius(X)
    return Ball2(c, r)
end

"""
    underapproximate(P::AbstractPolyhedron{N}, ::Type{Ellipsoid};
                     backend=default_sdp_solver(),
                     interior_point::AbstractVector{N}=chebyshev_center_radius(P)[1]
                    ) where {N<:Real}

Underapproximate a polyhedral set by the maximum-volume ellipsoid.

### Input

- `P`         -- polyhedral set
- `Ellipsoid` -- type for dispatch
- `backend`   -- (optional, default: `default_sdp_solver()`) backend to solve
                 the semidefinite program
- `interior_point` -- (optional, default: `chebyshev_center_radius(P)[1]`) an
                      interior point of the ellipsoid (needed for the algorithm:
                      see below for details)

### Output

An ellipsoid ``E`` such that ``E ⊆ P``. ``E`` has maximal volume if `P` is
bounded (no such guarantee is given in directions where `P` is unbounded).

### Notes

The maximum-volume ellipsoid is also called *Löwner-John ellipsoid*.

The algorithm requires to specify a point in the interior of the ellipsoid which
can be specified with the argument `interior_point`. If no such point is known
(i.e., `interior_point == nothing`), the implementation will compute the
Chebyshev center (see [`chebyshev_center_radius`](@ref)).

Note that the Chebyshev center will be computed using the default LP solver
and not the passed SDP solver `backend` (because it is not guaranteed that
`backend` supports LPs). If a custom LP solver should be used, the Chebyshev
center needs to be computed and passed manually.

### Algorithm

We use the package [`SetProg.jl`](https://github.com/blegat/SetProg.jl/) to
encode the problem directly.

An algorithm is described
[here](https://systemanalysisdpt-cmc-msu.github.io/ellipsoids/doc/chap_ellcalc.html#maximum-volume-ellipsoids).
"""
function underapproximate(P::AbstractPolyhedron{N}, ::Type{Ellipsoid};
                          backend=default_sdp_solver(),
                          interior_point::AbstractVector{N}=chebyshev_center_radius(P)[1]) where {N}
    require(@__MODULE__, :SetProg; fun_name="underapproximate")
    return _underapproximate_ellipsoid(P, Ellipsoid; backend=backend,
                                       interior_point=interior_point)
end

function load_setprog_underapproximate()
    return quote
        function _underapproximate_ellipsoid(P::AbstractPolyhedron{N}, ::Type{Ellipsoid};
                                             backend=default_sdp_solver(),
                                             interior_point::AbstractVector{N}=chebyshev_center_radius(P)[1]) where {N}
            # create SDP model
            model = Model(backend)

            # create symbolic ellipsoid S
            point = SetProg.InteriorPoint(interior_point)
            @variable(model, S, SetProg.Ellipsoid(point=point))

            # add subset constraint
            Q = polyhedron(P)  # convert P to a Polyhedra.polyhedron
            @constraint(model, S ⊆ Q)

            # add maximum-volume objective
            @objective(model, Max, SetProg.nth_root(SetProg.volume(S)))

            # solve SDP
            optimize!(model)

            # obtain solution
            set = SetProg.value(S)

            # convert to SetProg representation
            E = SetProg.Sets.ellipsoid(set)

            # convert result to a LazySets.Ellipsoid
            return convert(Ellipsoid, E)
        end

        function convert(::Type{<:Ellipsoid},
                         TE::SetProg.Sets.Translation{<:SetProg.Sets.Ellipsoid})
            return Ellipsoid(TE.c, inv(TE.set.Q))
        end

        function convert(::Type{<:Ellipsoid}, E::SetProg.Sets.Ellipsoid)
            return Ellipsoid(inv(E.Q))
        end
    end # quote
end # load_setprog_underapproximate
