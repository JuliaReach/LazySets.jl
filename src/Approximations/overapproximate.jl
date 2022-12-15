using LazySets: block_to_dimension_indices,
                substitute_blocks,
                fast_interval_pow,
                get_constrained_lowdimset

"""
    overapproximate(X::S, ::Type{S}, args...) where {S<:LazySet}

Overapproximating a set of type `S` with type `S` is a no-op.

### Input

- `X`       -- set
- `Type{S}` -- target set type
- `args`    -- further arguments (ignored)
- `kwargs`  -- further keyword arguments (ignored)

### Output

The input set.
"""
function overapproximate(X::S, ::Type{S}, args...; kwargs...) where {S<:LazySet}
    return X
end

"""
    overapproximate(S::LazySet)

Alias for `overapproximate(S, Hyperrectangle)` resp. `box_approximation(S)`.
"""
overapproximate(S::LazySet) = box_approximation(S)

"""
    overapproximate(S::LazySet, ::Type{<:Hyperrectangle})

Alias for `box_approximation(S)`.
"""
function overapproximate(S::LazySet, ::Type{<:Hyperrectangle})
    return box_approximation(S)
end

"""
    overapproximate(S::LazySet, ::Type{<:BallInf})

Alias for `ballinf_approximation(S)`.
"""
function overapproximate(S::LazySet, ::Type{<:BallInf})
    return ballinf_approximation(S)
end

"""
    overapproximate(S::LazySet{N},
                    ::Type{<:HPolygon},
                    [ε]::Real=Inf) where {N}

Overapproximate a given 2D set using iterative refinement.

### Input

- `S`        -- two-dimensional bounded set
- `HPolygon` -- target set type
- `ε`        -- (optional, default: `Inf`) error tolerance
- `prune`    -- (optional, default: `true`) flag for removing redundant
                constraints in the end

### Output

A polygon in constraint representation.

### Notes

The result is always a convex overapproximation of the input set.

If no error tolerance ε is given, or is `Inf`, the result is a box-shaped
polygon. For convex input sets, the result is an ε-close polygonal
overapproximation with respect to the Hausdorff distance.
"""
function overapproximate(S::LazySet{N},
                         ::Type{<:HPolygon},
                         ε::Real=Inf;
                         prune::Bool=true) where {N}
    @assert dim(S) == 2 "epsilon-close overapproximation is only available " *
                        "for two-dimensional sets"
    if ε == Inf
        constraints = Vector{HalfSpace{N, Vector{N}}}(undef, 4)
        constraints[1] = HalfSpace(DIR_EAST(N), ρ(DIR_EAST(N), S))
        constraints[2] = HalfSpace(DIR_NORTH(N), ρ(DIR_NORTH(N), S))
        constraints[3] = HalfSpace(DIR_WEST(N), ρ(DIR_WEST(N), S))
        constraints[4] = HalfSpace(DIR_SOUTH(N), ρ(DIR_SOUTH(N), S))
        return HPolygon(constraints, sort_constraints=false)
    else
        P = overapproximate_hausdorff(S, ε)
        if prune
            remove_redundant_constraints!(P)
        end
        return P
    end
end

# error messages for dispatch on HPolytope
function overapproximate(X::S, ::Type{<:HPolytope}) where {S<:LazySet}
    if dim(X) == 2
        throw(ArgumentError("epsilon-close approximation is only available " *
              "using polygons in constraint representation; try " *
              "`overapproximate(X, HPolygon)`"))
    else
        throw(ArgumentError("epsilon-close approximation is only available " *
              "for two-dimensional sets; try " *
              "`overapproximate(X, HPolytope, dirs)` where `dirs` are " *
              "template directions, e.g., `BoxDirections` or `OctDirections`"))
    end
end

"""
    overapproximate(S::LazySet, ε::Real)

Alias for `overapproximate(S, HPolygon, ε)`.
"""
function overapproximate(S::LazySet, ε::Real)
    return overapproximate(S, HPolygon, ε)
end

# special case: overapproximation of empty set
overapproximate(∅::EmptySet, args...; kwargs...) = ∅

# disambiguation
for ST in LazySets.subtypes(LazySet, true)
    if ST == HPolygon  # must be defined separately below with extra argument
        continue
    end
    @eval overapproximate(∅::EmptySet, ::Type{<:$ST}) = ∅
end
overapproximate(∅::EmptySet, ::Real) = ∅
overapproximate(∅::EmptySet, ::Type{<:HPolygon}, ::Real=Inf; kwargs...) = ∅
overapproximate(∅::EmptySet, ::Type{<:EmptySet}, args...; kwargs...) = ∅

"""
    overapproximate(X::LazySet{N}, dir::AbstractDirections;
                    [prune]::Bool=true) where {N}

Overapproximate a (possibly unbounded) set with template directions.

### Input

- `X`     -- set
- `dir`   -- directions
- `prune` -- (optional, default: `true`) flag for removing redundant constraints

### Output

A polyhedron overapproximating the set `X` with the directions from `dir`. The
overapproximation is computed using the support function.
The result is an `HPolytope` if it is bounded and otherwise an `HPolyhedron`.
"""
function overapproximate(X::LazySet, dir::AbstractDirections; prune::Bool=true)
    H = _overapproximate_directions(X, dir, prune)

    # if input is bounded and directions are bounding => output is bounded
    # otherwise, check boundedness of the output
    if (isbounded(X) && isbounding(dir)) || _isbounded_stiemke(H)
        return HPolytope(H, check_boundedness=false)
    else
        return HPolyhedron(H)
    end
end

function _overapproximate_directions(X::LazySet{N},
                                     dir::AbstractDirections{N, VN},
                                     prune::Bool) where {N, VN}
    H = Vector{HalfSpace{N, VN}}()
    sizehint!(H, length(dir))

    for d in dir
        sf = ρ(d, X)
        if !isinf(sf)
            push!(H, HalfSpace(d, sf))
        end
    end
    if prune
        remove_redundant_constraints!(H) || throw(ArgumentError("unable to " *
                                                "remove redundant constraints"))
    end
    return H
end

# alias with HPolytope type as second argument
function overapproximate(X::LazySet, ::Type{<:HPolytope},
                         dirs::AbstractDirections; prune::Bool=true)
    P = overapproximate(X, dirs, prune=prune)
    P isa HPolytope || throw(ArgumentError("cannot overapproximate with an " *
        "`HPolytope` because the set is unbounded; try using an `HPolyhedron`"))
    return P
end

# alias with HPolyhedron type as second argument
function overapproximate(X::LazySet, ::Type{<:HPolyhedron},
                         dirs::AbstractDirections; prune::Bool=true)
    H = _overapproximate_directions(X, dir, prune)
    return HPolyhedron(H)
end

# disambiguation
overapproximate(∅::EmptySet, ::Type{<:HPolytope}, dirs::AbstractDirections;
                prune::Bool=true) = ∅
overapproximate(∅::EmptySet, ::Type{<:HPolyhedron}, dirs::AbstractDirections;
                prune::Bool=true) = ∅

"""
    overapproximate(X::LazySet{N}, dir::Type{<:AbstractDirections}) where {N}

Overapproximate a set with template directions.

### Input

- `X`   -- set
- `dir` -- type of direction representation

### Output

A polyhedron overapproximating the set `X` with the directions from `dir`.
The result is an `HPolytope` if it is bounded and otherwise an `HPolyhedron`.
"""
function overapproximate(X::LazySet{N},
                         dir::Type{<:AbstractDirections}; kwargs...) where {N}
    return overapproximate(X, dir{N}(dim(X)); kwargs...)
end

# disambiguation
overapproximate(∅::EmptySet, dir::Type{<:AbstractDirections}; kwargs...) = ∅
overapproximate(∅::EmptySet, dir::AbstractDirections; prune::Bool=true) = ∅

function overapproximate_cap_helper(X::LazySet,
                                    P::AbstractPolyhedron,  # polyhedron
                                    dir::AbstractDirections;
                                    kwargs...
                                   )
    Hi = constraints_list(P)
    m = length(Hi)
    N = promote_type(eltype(X), eltype(P))
    constraints = Vector{HalfSpace{N, Vector{N}}}() # TODO: use directions type, see #2031
    sizehint!(constraints, length(dir))
    return_type = HPolytope

    for di in dir
        ρ_X_Hi_min = ρ(di, X ∩ Hi[1], kwargs...)
        for i in 2:m
            ρ_X_Hi = ρ(di, X ∩ Hi[i], kwargs...)
            if ρ_X_Hi < ρ_X_Hi_min
                ρ_X_Hi_min = ρ_X_Hi
            end
        end
        if isinf(ρ_X_Hi_min)
            # unbounded in this direction => return a polyhedron later
            return_type = HPolyhedron
        else
            push!(constraints, _normal_Vector(HalfSpace(di, ρ_X_Hi_min))) # TODO: remove vector
        end
    end
    return return_type(constraints)
end

"""
    overapproximate(cap::Intersection{N, <:LazySet, <:AbstractPolyhedron},
                    dir::AbstractDirections;
                    kwargs...
                   ) where {N}

Overapproximate the intersection between a set and a polyhedron given a set of
template directions.

### Input

- `cap`         -- intersection of a set and a polyhedron
- `dir`         -- template directions
- `kwargs`      -- additional arguments that are passed to the support function
                   algorithm

### Output

A polytope or polyhedron in H-representation such that the normal direction of
each half-space is given by an element of `dir`.

### Algorithm

Let `di` be a direction drawn from the set of template directions `dir`.
Let `X` be the set and let `P` be the polyhedron. We overapproximate the set
`X ∩ H` with a polytope or polyhedron in constraint representation using a given
set of template directions `dir`.

The idea is to solve the univariate optimization problem `ρ(di, X ∩ Hi)` for
each half-space of the set `P` and then take the minimum.
This gives an overapproximation of the exact support function.

This algorithm is inspired from [1].

### Notes

This method relies on having available the `constraints_list` of the polyhedron
`P`.

This method may return a non-empty set even if the original set is empty.

[1] G. Frehse, R. Ray. *Flowpipe-Guard Intersection for Reachability
Computations with Support Functions*. ADHS 2012.
"""
function overapproximate(cap::Intersection{N,  # TODO use better mechanism to detect polyhedral set
                                           <:LazySet,
                                           <:AbstractPolyhedron},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron,
                                           <:LazySet},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.Y, cap.X, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron,
                                           <:AbstractPolyhedron},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    # important: the result may not be a polytope!
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope,
                                           <:AbstractPolyhedron},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron,
                                           <:AbstractPolytope},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.Y, cap.X, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope,
                                           <:AbstractPolytope},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

"""
    overapproximate(cap::Intersection{N, <:HalfSpace, <:AbstractPolytope},
                    dir::AbstractDirections;
                    [kwargs]...
                   ) where {N}

Overapproximate the intersection between a half-space and a polytope given a set
of template directions.

### Input

- `cap`         -- intersection of a half-space and a polytope
- `dir`         -- template directions
- `kwargs`      -- additional arguments that are passed to the support function
                   algorithm

### Output

A polytope in H-representation such that the normal direction of each half-space
is given by an element of `dir`.
"""
function overapproximate(cap::Intersection{N,
                                           <:HalfSpace,
                                           <:AbstractPolytope},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    H = HPolytope{N, Vector{N}}()
    c = H.constraints
    push!(c, _normal_Vector(cap.X))
    append!(c, _normal_Vector(cap.Y))
    return overapproximate(H, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope,
                                           <:HalfSpace},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate(swap(cap), dir; kwargs...)
end

function overapproximate(P::SimpleSparsePolynomialZonotope,
                         ::Type{<:UnionSetArray{Zonotope}};
                         nsdiv=10, partition=nothing)
    q = nparams(P)
    dom = IA.IntervalBox(IA.Interval(-1, 1), q)
    cells = IA.mince(dom, isnothing(partition) ? nsdiv : partition)
    return UnionSetArray([overapproximate(P, Zonotope, c) for c in cells])
end

"""
    overapproximate(Z::AbstractZonotope, ::Type{<:HParallelotope},
                    [indices]=1:dim(Z))

Overapproximate a zonotopic set with a parallelotope in constraint
representation.

### Input

- `Z`              -- zonotopic set
- `HParallelotope` -- target set type
- `indices`        -- (optional; default: `1:dim(Z)`) generator indices selected
                       when constructing the parallelotope

### Output

An overapproximation of the given zonotopic set using a parallelotope.

### Algorithm

The algorithm is based on Proposition 8 discussed in Section 5 of [1].

[1] Althoff, M., Stursberg, O., & Buss, M. (2010). *Computing reachable sets of
hybrid systems using a combination of zonotopes and polytopes*. Nonlinear
analysis: hybrid systems, 4(2), 233-249.
"""
function overapproximate(Z::AbstractZonotope, ::Type{<:HParallelotope},
                         indices=1:dim(Z))
    Zred = _overapproximate_hparallelotope(Z, indices)
    return convert(HParallelotope, Zred)
end

function _overapproximate_hparallelotope(Z::AbstractZonotope, indices=1:dim(Z))
    length(indices) == dim(Z) || throw(ArgumentError("the number of " *
        "generator indices is $(length(indices)), but it was expected to be " *
        "$(dim(Z))"))

    p, n = ngens(Z), dim(Z)
    if p == n
        return Z
    elseif p < n
        error("the zonotope order is $(order(Z)) but it should be at least 1")
    end

    G = genmat(Z)
    Γ = G[:, indices]
    □Γ⁻¹Z = box_approximation(linear_map(inv(Γ), Z))
    return linear_map(Γ, □Γ⁻¹Z)
end

"""
    overapproximate(X::Intersection{N, <:AbstractZonotope, <:Hyperplane},
                    dirs::AbstractDirections) where {N}

Overapproximate the intersection between a zonotopic set and a hyperplane with a
polyhedron or polytope using the given directions.

### Input

- `X`    -- intersection between a zonotopic set and a hyperplane
- `dirs` -- type of direction representation

### Output

An overapproximation of the intersection between a zonotopic set and a
hyperplane. If the directions are bounding, the result is an `HPolytope`,
otherwise the result is an `HPolyhedron`.

### Algorithm

This function implements [Algorithm 8.1, 1].

[1] Colas Le Guernic. *Reachability Analysis of Hybrid Systems with Linear
continuous dynamics* (Doctoral dissertation). 2009.
"""
function overapproximate(X::Intersection{N, <:AbstractZonotope, <:Hyperplane},
                         dirs::AbstractDirections) where {N}
    dim(X) == dim(dirs) || throw(ArgumentError("the dimension of the set " *
        "$(dim(X)) does not match the dimension of the directions $(dim(dirs))"))
    Z, G = X.X, X.Y

    if isdisjoint(Z, G)
        return EmptySet{N}(dim(Z))
    end

    n = G.a  # normal vector to the hyperplane
    γ = G.b  # displacement of the hyperplane
    Lᵧ = Line2D([one(N), zero(N)], γ)  # line (x, y) : x = γ

    constraints = Vector{HalfSpace{N, eltype(dirs)}}()
    for l in dirs
        Πₙₗ = vcat(n', l')  # projection map
        πZₙₗ = linear_map(Πₙₗ, Z)

        ρₗ = LazySets._bound_intersect_2D(πZₙₗ, Lᵧ)

        push!(constraints, HalfSpace(l, ρₗ))
    end
    T = isbounding(dirs) ? HPolytope : HPolyhedron
    return T(constraints)
end

# symmetric method
function overapproximate(X::Intersection{N, <:Hyperplane, <:AbstractZonotope},
                         dirs::AbstractDirections) where {N}
    return overapproximate(X.Y ∩ X.X, dirs)
end

# overload on direction type
function overapproximate(X::Intersection{N, <:AbstractZonotope,
                         <:Hyperplane}, dirs::Type{<:AbstractDirections}) where {N}
    return overapproximate(X, dirs(dim(X)))
end

# symmetric method
function overapproximate(X::Intersection{N, <:Hyperplane, <:AbstractZonotope},
                         dirs::Type{<:AbstractDirections}) where {N}
    return overapproximate(X.Y ∩ X.X, dirs(dim(X)))
end

"""
    overapproximate(QM::QuadraticMap{N, <:SparsePolynomialZonotope},
                    ::Type{<:SparsePolynomialZonotope}) where {N}

Overapproximate a quadratic map of a sparse polynomial zonotope with a sparse
polynomial zonotope.

### Input

- `QM`                       -- quadratic map of a sparse polynomial zonotope
- `SparsePolynomialZonotope` -- target type

### Output

A sparse polynomial zonotope overapproximating the quadratic map of a sparse
polynomial zonotope.

### Algorithm

This method implements Proposition 13 of [1].

[1] N. Kochdumper and M. Althoff. *Sparse Polynomial Zonotopes: A Novel Set
Representation for Reachability Analysis*. Transactions on Automatic Control
2021.
"""
function overapproximate(QM::QuadraticMap{N, <:SparsePolynomialZonotope},
                         ::Type{<:SparsePolynomialZonotope}) where {N}
    PZ = QM.X
    Q = QM.Q

    p = nparams(PZ)
    q = ngens_indep(PZ)

    c = center(PZ)
    G = genmat_dep(PZ)
    GI = genmat_indep(PZ)
    Ĝ = hcat(G, GI)
    Ê = cat(expmat(PZ), Matrix(1 * I, q, q); dims=(1, 2))

    PZ1 = SimpleSparsePolynomialZonotope(c, Ĝ, Ê)
    PZbar = quadratic_map(Q, PZ1)
    c̄ = center(PZbar)
    Ē = expmat(PZbar)
    Ḡ = genmat(PZbar)

    K = [!iszero(Ei[p+1:end]) for Ei in eachcol(Ē)]
    H = .!K
    PZK = SimpleSparsePolynomialZonotope(c̄, Ḡ[:, K], Ē[:, K])
    Z = overapproximate(PZK, Zonotope)
    cz = center(Z)
    Gz = genmat(Z)
    return SparsePolynomialZonotope(cz, Ḡ[:, H], Gz, Ē[1:p, H], indexvector(PZ))
end

# function to be loaded by Requires
function load_paving_overapproximation()
return quote

using .IntervalConstraintProgramming: Paving

"""
    overapproximate(p::Paving{L, N}, dirs::AbstractDirections{N, VN})
        where {L, N, VN}

Overapproximate a Paving-type set representation with a polyhedron in constraint
representation.

### Input

- `p`    -- paving
- `dirs` -- template directions

### Output

An overapproximation of a paving using a polyhedron in constraint representation
(`HPolyhedron`) with constraints in direction `dirs`.

### Algorithm

This function takes the union of the elements at the boundary of `p`, first
converted into hyperrectangles, and then calculates the support function of the
set along each  direction in dirs, to compute the `HPolyhedron` constraints.

### Requires IntervalConstraintProgramming
"""
function overapproximate(p::Paving{L, N}, dirs::AbstractDirections{N, VN}) where {L, N, VN}
    # enclose outer approximation
    Uouter = UnionSetArray(convert.(Hyperrectangle, p.boundary))
    constraints = [HalfSpace(d, ρ(d, Uouter)) for d in dirs]
    return HPolyhedron(constraints)
end

# alias with HPolyhedron type as second argument
function overapproximate(p::Paving{L, N}, ::Type{<:HPolyhedron},
                         dirs::AbstractDirections{N, VN}) where {L, N, VN}
    return overapproximate(p, dirs)
end

end end  # quote / load_paving_overapproximation
