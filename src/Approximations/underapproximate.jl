"""
    underapproximate(X::ConvexSet{N}, dirs::AbstractDirections;
                    [apply_convex_hull]::Bool=false) where {N}

Compute the underapproximation of a convex set by sampling support vectors.

### Input

- `X`                 -- set
- `dirs`              -- directions
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the support vectors with a convex hull operation
### Output

The `VPolytope` obtained by taking the convex hull of the support vectors of `X`
along the directions determined by `dirs`.

### Notes

Since the support vectors are not always unique, this algorithm may return
a strict underapproximation even if the set can be exactly approximated using
the given template.
"""
function underapproximate(X::ConvexSet{N}, dirs::AbstractDirections;
                          apply_convex_hull::Bool=false) where {N}
    vinner = Vector{Vector{N}}(undef, length(dirs))
    underapproximate!(vinner, X, dirs; apply_convex_hull=apply_convex_hull)
end

# in-place version
function underapproximate!(vinner, X::ConvexSet{N}, dirs::AbstractDirections;
                           apply_convex_hull::Bool=false) where {N}
   @assert dim(X) == dim(dirs) "the dimension of the set, $(dim(X)), doesn't match " *
                               "the dimension of the template directions, $(dim(dirs))"

    @inbounds for (i, di) in enumerate(dirs)
        vinner[i] = σ(di, X)
    end
    if apply_convex_hull
        convex_hull!(vinner)
    end
    return VPolytope(vinner)
end

function underapproximate(X::ConvexSet, dirs::Type{<:AbstractDirections}; kwargs...)
    return underapproximate(X, dirs(dim(X)), kwargs...)
end

"""
    underapproximate(X::ConvexSet, ::Type{<:Hyperrectangle};
                     solver=nothing) where {N}

Underapproximate a polygon with a hyperrectangle of maximal area.

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
function underapproximate(X::ConvexSet{N}, ::Type{<:Hyperrectangle};
                          solver=nothing) where {N}
    require(:Ipopt; fun_name="underapproximate")

    solver = default_nln_solver(N)
    return _underapproximate_box(X, solver)
end
