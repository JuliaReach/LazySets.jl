"""
    kamenev(S::LazySet, tol::Real; [maxiter]::Int=100)

Compute an inner and outer polyhedral approximation of a convex set.

### Input

- `S`       -- convex set
- `tol`     -- error tolerance (Hausdorff distance between inner and outer approximation)
- `maxiter` -- (optional, default: `100`) maximum number of iterations

### Output

A tuple `(inner, outer, error)` where `inner` is an inner approximation, `outer` is an outer approximation,
and `error` is the final error.
"""
function kamenev(S::LazySet, tol; maxiter = 100)
    # v representation of an inner approximation
    inner_V = underapproximate(S,BoxDirections(dim(S)))
    inner_H = convert(HPolytope,inner_V) # H representation of the inner polytope

    # in this initial loop, we will construct the outer approximation
    spaces = [HalfSpace(face.a, ρ(face.a,S)) for face in constraints_list(inner_H)]
    outer = HPolytope(spaces) # the outer approximation

    cur_err = Inf
    for iter in 1:maxiter
        (cur_err,add_ind) = findmax(constraints_list(inner_H)) do face
            ρ(face.a,outer)-face.b
        end
        
        cur_err < tol && break
        worst_face = constraints_list(inner_H)[add_ind]
        
        push!(inner_V.vertices,σ(worst_face.a,S)) # this point must be added to the list of inner vertices
        inner_V = VPolytope(inner_V.vertices)
        inner_H = convert(HPolytope,inner_V) # this isn't ideal. inner_H was still largely correct

        outer = concretize(outer ∩ HalfSpace(worst_face.a,ρ(worst_face.a,S)))
    end

    return (inner_V, outer, cur_err)
end

"""
    underapproximate(S::LazySet, ::Type{<:VPolytope}, tol::Real; kwargs...)

Underapproximate a convex set with a `VPolytope` using Kamenev's algorithm.

### Input

- `S`    -- convex set
- `type` -- target type `VPolytope`
- `tol`  -- error tolerance

### Output

A `VPolytope` underapproximation.
"""
underapproximate(S::LazySet,::Type{<:VPolytope}, args...; kwargs...) = kamenev(S,args...; kwargs...)[1]

"""
    overapproximate(S::LazySet, ::Type{<:HPolytope}, tol::Real; kwargs...)

Overapproximate a convex set with an `HPolytope` using Kamenev's algorithm.

### Input

- `S`    -- convex set
- `type` -- target type `HPolytope`
- `tol`  -- error tolerance

### Output

An `HPolytope` overapproximation.
"""
overapproximate(S::LazySet,::Type{<:HPolytope}, args...; kwargs...) = kamenev(S,args...; kwargs...)[2]
