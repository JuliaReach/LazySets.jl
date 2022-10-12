"""
    underapproximate(X::S, dirs::AbstractDirections;
                    [apply_convex_hull]::Bool=false) where {N, S<:LazySet{N}}

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
function underapproximate(X::S, dirs::AbstractDirections;
                          apply_convex_hull::Bool=false) where {N, S<:LazySet{N}}
    if !isconvextype(S)
        error("this underapproximation is only available for convex sets")
    end
    @assert dim(X) == dim(dirs) "the dimension of the set, $(dim(X)), does " *
              "not match the dimension of the template directions, $(dim(dirs))"

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

The algorithm is taken from [1, Theorem 17] and solves a convex program (in fact
a linear program with nonlinear objective). (The objective is modified to an
equivalent version due to solver issues.)

[1] Mehdi Behroozi - *Largest inscribed rectangles in geometric convex sets.*
arXiv:1905.13246.
"""
function underapproximate(X::LazySet{N}, ::Type{<:Hyperrectangle};
                          solver=nothing) where {N}
    require(@__MODULE__, :Ipopt; fun_name="underapproximate")

    if isnothing(solver)
        solver = default_nln_solver(N)
    end
    return _underapproximate_box(X, solver)
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
