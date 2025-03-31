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
    overapproximate(S::LazySet, T::Type{<:LazySet}, [args]...; [kwargs]...)

Default overapproximation method that falls back to `convert`.

### Input

- `X`       -- set
- `Type{S}` -- target set type
- `args`    -- further arguments
- `kwargs`  -- further keyword arguments

### Output

The result of `convert`, or a `MethodError` if no such method exists.
"""
function overapproximate(X::T1, T2::Type{<:LazySet}, args...; kwargs...) where {T1<:LazySet}
    return convert(T2, X, args...; kwargs...)
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
        constraints = Vector{HalfSpace{N,Vector{N}}}(undef, 4)
        constraints[1] = HalfSpace(DIR_EAST(N), ρ(DIR_EAST(N), S))
        constraints[2] = HalfSpace(DIR_NORTH(N), ρ(DIR_NORTH(N), S))
        constraints[3] = HalfSpace(DIR_WEST(N), ρ(DIR_WEST(N), S))
        constraints[4] = HalfSpace(DIR_SOUTH(N), ρ(DIR_SOUTH(N), S))
        return HPolygon(constraints; sort_constraints=false)
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
    if dim(S) == 1
        return overapproximate(S, Interval)
    else
        return overapproximate(S, HPolygon, ε)
    end
end

# special case: overapproximation of empty set
overapproximate(∅::EmptySet, args...; kwargs...) = ∅

"""
    overapproximate(X::LazySet{N}, dirs::AbstractDirections;
                    [prune]::Bool=true) where {N}

Overapproximate a (possibly unbounded) set with template directions.

### Input

- `X`     -- set
- `dirs`  -- directions
- `prune` -- (optional, default: `true`) flag for removing redundant constraints

### Output

A polyhedron overapproximating the set `X` with the directions from `dirs`. The
overapproximation is computed using the support function.
The result is an `HPolytope` if it is bounded and otherwise an `HPolyhedron`.
"""
function overapproximate(X::LazySet, dirs::AbstractDirections; prune::Bool=true)
    H = _overapproximate_directions(X, dirs, prune)

    # if input is bounded and directions are bounding => output is bounded
    # otherwise, check boundedness of the output
    if (isbounded(X) && isbounding(dirs)) || _isbounded_stiemke(H)
        return HPolytope(H; check_boundedness=false)
    else
        return HPolyhedron(H)
    end
end

function _overapproximate_directions(X::LazySet{N},
                                     dirs::AbstractDirections{N,VN},
                                     prune::Bool) where {N,VN}
    H = Vector{HalfSpace{N,VN}}()
    sizehint!(H, length(dirs))

    for d in dirs
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
    P = overapproximate(X, dirs; prune=prune)
    P isa HPolytope || throw(ArgumentError("cannot overapproximate with an " *
                                           "`HPolytope` because the set is unbounded; try using an `HPolyhedron`"))
    return P
end

# alias with HPolyhedron type as second argument
function overapproximate(X::LazySet, ::Type{<:HPolyhedron},
                         dirs::AbstractDirections; prune::Bool=true)
    H = _overapproximate_directions(X, dirs, prune)
    return HPolyhedron(H)
end

"""
    overapproximate(X::LazySet{N}, dirs::Type{<:AbstractDirections}) where {N}

Overapproximate a set with template directions.

### Input

- `X`    -- set
- `dirs` -- type of direction representation

### Output

A polyhedron overapproximating the set `X` with the directions from `dirs`.
The result is an `HPolytope` if it is bounded and otherwise an `HPolyhedron`.
"""
function overapproximate(X::LazySet{N},
                         dirs::Type{<:AbstractDirections}; kwargs...) where {N}
    return overapproximate(X, dirs{N}(dim(X)); kwargs...)
end

function overapproximate_cap_helper(X::LazySet,
                                    P::AbstractPolyhedron,  # polyhedron
                                    dirs::AbstractDirections;
                                    kwargs...)
    Hi = constraints_list(P)
    m = length(Hi)
    N = promote_type(eltype(X), eltype(P))
    constraints = Vector{HalfSpace{N,Vector{N}}}() # TODO: use directions type, see #2031
    sizehint!(constraints, length(dirs))
    return_type = HPolytope

    for di in dirs
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
                    dirs::AbstractDirections;
                    kwargs...
                   ) where {N}

Overapproximate the intersection between a set and a polyhedron given a set of
template directions.

### Input

- `cap`    -- intersection of a set and a polyhedron
- `dirs`   -- template directions
- `kwargs` -- additional arguments that are passed to the support function
              algorithm

### Output

A polytope or polyhedron in H-representation such that the normal direction of
each half-space is given by an element of `dirs`.

### Algorithm

Let `di` be a direction drawn from the set of template directions `dirs`.
Let `X` be the set and let `P` be the polyhedron. We overapproximate the set
`X ∩ H` with a polytope or polyhedron in constraint representation using a given
set of template directions `dirs`.

The idea is to solve the univariate optimization problem `ρ(di, X ∩ Hi)` for
each half-space of the set `P` and then take the minimum.
This gives an overapproximation of the exact support function.

This algorithm is inspired from [Frehse012](@citet).

### Notes

This method relies on having available the `constraints_list` of the polyhedron
`P`.

This method may return a non-empty set even if the original set is empty.
"""
function overapproximate(cap::Intersection{N,  # TODO use better mechanism to detect polyhedral set
                                           <:LazySet,
                                           <:AbstractPolyhedron},
                         dirs::AbstractDirections;
                         kwargs...) where {N}
    return overapproximate_cap_helper(first(cap), second(cap), dirs; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron,
                                           <:LazySet},
                         dirs::AbstractDirections;
                         kwargs...) where {N}
    return overapproximate_cap_helper(second(cap), first(cap), dirs; kwargs...)
end

"""
    overapproximate(cap::Intersection{N, <:HalfSpace, <:AbstractPolytope},
                    dirs::AbstractDirections;
                    [kwargs]...
                   ) where {N}

Overapproximate the intersection between a half-space and a polytope given a set
of template directions.

### Input

- `cap`    -- intersection of a half-space and a polytope
- `dirs`   -- template directions
- `kwargs` -- additional arguments that are passed to the support function
              algorithm

### Output

A polytope in H-representation such that the normal direction of each half-space
is given by an element of `dirs`.
"""
function overapproximate(cap::Intersection{N,
                                           <:HalfSpace,
                                           <:AbstractPolytope},
                         dirs::AbstractDirections;
                         kwargs...) where {N}
    H = HPolytope{N,Vector{N}}()
    c = H.constraints
    push!(c, _normal_Vector(first(cap)))
    append!(c, _normal_Vector(second(cap)))
    return overapproximate(H, dirs; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope,
                                           <:HalfSpace},
                         dirs::AbstractDirections;
                         kwargs...) where {N}
    return overapproximate(swap(cap), dirs; kwargs...)
end

function overapproximate(P::SimpleSparsePolynomialZonotope,
                         ::Type{<:UnionSetArray{Zonotope}};
                         nsdiv=10, partition=nothing)
    q = nparams(P)
    dom = IA.IntervalBox(IA.interval(-1, 1), q)
    cells = IA.mince(dom, isnothing(partition) ? nsdiv : partition)
    return UnionSetArray([overapproximate(P, Zonotope, c) for c in cells])
end

function overapproximate(P::SparsePolynomialZonotope,
                         T::Type{<:UnionSetArray{Zonotope}};
                         nsdiv=10, partition=nothing)
    SSPZ = convert(SimpleSparsePolynomialZonotope, P)
    return overapproximate(SSPZ, T; nsdiv=nsdiv, partition=partition)
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

The algorithm is based on Proposition 8 discussed in [AlthoffSB10; Section 5](@citet).
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

This function implements [LeGuernic09; Algorithm 8.1](@citet).
"""
function overapproximate(X::Intersection{N,<:AbstractZonotope,<:Hyperplane},
                         dirs::AbstractDirections) where {N}
    dim(X) == dim(dirs) || throw(ArgumentError("the dimension of the set " *
                                               "$(dim(X)) does not match the dimension of the directions $(dim(dirs))"))
    Z, G = first(X), second(X)

    if isdisjoint(Z, G)
        return EmptySet{N}(dim(Z))
    end

    n = G.a  # normal vector to the hyperplane
    γ = G.b  # displacement of the hyperplane
    Lᵧ = Line2D([one(N), zero(N)], γ)  # line (x, y) : x = γ

    constraints = Vector{HalfSpace{N,eltype(dirs)}}()
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
function overapproximate(X::Intersection{N,<:Hyperplane,<:AbstractZonotope},
                         dirs::AbstractDirections) where {N}
    return overapproximate(second(X) ∩ first(X), dirs)
end

# overload on direction type
function overapproximate(X::Intersection{N,<:AbstractZonotope,
                                         <:Hyperplane}, dirs::Type{<:AbstractDirections}) where {N}
    return overapproximate(X, dirs(dim(X)))
end

# symmetric method
function overapproximate(X::Intersection{N,<:Hyperplane,<:AbstractZonotope},
                         dirs::Type{<:AbstractDirections}) where {N}
    return overapproximate(second(X) ∩ first(X), dirs(dim(X)))
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

This method implements [KochdumperA21; Proposition 13](@citet).
"""
function overapproximate(QM::QuadraticMap{N,<:SparsePolynomialZonotope},
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

    K = [!iszero(Ei[(p + 1):end]) for Ei in eachcol(Ē)]
    H = .!K
    PZK = SimpleSparsePolynomialZonotope(c̄, Ḡ[:, K], Ē[:, K])
    Z = overapproximate(PZK, Zonotope)
    cz = center(Z)
    Gz = genmat(Z)
    return SparsePolynomialZonotope(cz, Ḡ[:, H], Gz, Ē[1:p, H], indexvector(PZ))
end

"""
    overapproximate(P::VPolygon, ::Type{<:LinearMap{N,<:Hyperrectangle}}) where {N}

Overapproximate a convex polygon with a minimal-area rotated rectangle.

### Input

- `P` -- convex polygon
- `LinearMap{N,<:Hyperrectangle}` -- target type

### Output

A `LinearMap` of a `Hyperrectangle`.

### Algorithm

This method follows an approach described in
[this post](https://gis.stackexchange.com/a/22934), which itself is based on
[this post](https://gis.stackexchange.com/a/22904).

Generally, the idea is that the rotated rectangle must share at least one edge
with the polygon. Thus, it suffices to try out finitely many rectangles. Some
tricks from linear algebra allow to construct the corresponding rotations and
rectangles elegantly.
"""
function overapproximate(P::VPolygon, ::Type{<:LinearMap{N,<:Hyperrectangle}}) where {N}
    min_area = N(Inf)
    vert_P = P.vertices
    m = length(vert_P)

    center = Vector{N}(undef, 2)
    radius = Vector{N}(undef, 2)
    R = Matrix{N}(undef, 2, 2)

    # iterate over all polygon edges
    @inbounds for i in eachindex(vert_P)
        # edge (v_i, v_{i+1})
        a = vert_P[i]
        next_idx = i % m + 1
        b = vert_P[next_idx]
        e = b - a
        v = normalize(e)
        w = [-v[2], v[1]]

        min_x, max_x = extrema(dot(vertex, v) for vertex in vert_P)
        min_y, max_y = extrema(dot(vertex, w) for vertex in vert_P)

        current_area = (max_x - min_x) * (max_y - min_y)

        if current_area < min_area
            min_area = current_area
            center .= ((max_x + min_x) / 2, (max_y + min_y) / 2)
            radius .= ((max_x - min_x) / 2, (max_y - min_y) / 2)
            R[:, 1] = v
            R[:, 2] = w
        end
    end

    min_rectangle = Hyperrectangle(center, radius)

    return LinearMap(R, min_rectangle)
end

"""
    overapproximate(P::SparsePolynomialZonotope{N}, ::Type{<:VPolytope}) where {N}


Overapproximate a sparse polynomial zonotope with a polytope in vertex representation.

### Input

- `P`         -- sparse polynomial zonotope
- `VPolytope` -- target type

### Output

A `VPolytope` that overapproximates the sparse polynomial zonotope.

### Algorithm

This method implements [Kochdumper21a; Proposition 3.1.15](@citet).
The idea is to split `P` into a linear and nonlinear part (such that `P = P₁ ⊕ P₂`). 
The nonlinear part is enclosed by a zonotope. Then we combine the vertices
of both sets and finally apply a convex-hull algorithm.
"""
function overapproximate(P::SparsePolynomialZonotope{N}, ::Type{<:VPolytope}) where {N}
    c = center(P)
    G = genmat_dep(P)
    GI = genmat_indep(P)
    E = expmat(P)
    idx = P.idx

    H = [j for j in 1:size(E, 2) if any(E[:, j] .> 1)]
    K = setdiff(1:size(E, 2), H)

    if !isempty(H)
        SPZ₂ = SparsePolynomialZonotope(c, G[:, H], zeros(N, length(c), 0), E[:, H], idx)
        Z = overapproximate(SPZ₂, Zonotope)
        c_z = center(Z)
        GI_mod = hcat(GI, genmat(Z))
    else
        c_z = c
        GI_mod = GI
    end

    G_mod = G[:, K]
    E_mod = E[:, K]
    
    # P̄ = SparsePolynomialZonotope(c_z, G_mod, GI_mod, E_mod, idx)
    # Compute vertices of a Z-representation
    p = size(E, 1)
    dep_params = Iterators.product(fill([-one(N), one(N)], p)...)
    indep_params = Iterators.product(fill([-one(N), one(N)], size(GI_mod, 2))...)

    V = Vector{Vector{N}}()
    for α in dep_params
        dep_term = zeros(N, size(c))
        for j in axes(G_mod, 2)
            prod = one(N)
            for k in 1:p
                if E_mod[k, j] == 1
                    prod *= α[k]
                end
            end
            dep_term += prod * G_mod[:, j]
        end
        for β in indep_params
            indep_term = zeros(N, size(c))
            for j in axes(GI_mod, 2)
                indep_term += β[j] * GI_mod[:, j]
            end
            point = c_z + dep_term + indep_term
            push!(V, point)
        end
    end

    convex_hull!(V)

    return VPolytope(V)
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

        This algorithm requires the IntervalConstraintProgramming package.
        """
        function overapproximate(p::Paving{L,N}, dirs::AbstractDirections{N,VN}) where {L,N,VN}
            # enclose outer approximation
            Uouter = UnionSetArray(convert.(Hyperrectangle, p.boundary))
            constraints = [HalfSpace(d, ρ(d, Uouter)) for d in dirs]
            return HPolyhedron(constraints)
        end

        # alias with HPolyhedron type as second argument
        function overapproximate(p::Paving{L,N}, ::Type{<:HPolyhedron},
                                 dirs::AbstractDirections{N,VN}) where {L,N,VN}
            return overapproximate(p, dirs)
        end
    end
end  # quote / load_paving_overapproximation

# ============== #
# disambiguation #
# ============== #

for ST in subtypes(LazySet, true)
    if ST == HPolygon  # must be defined separately below with extra argument
        continue
    end
    @eval overapproximate(∅::EmptySet, ::Type{<:$ST}) = ∅
end

overapproximate(∅::EmptySet, ::Type{<:LazySet}, args...; kwargs...) = ∅

overapproximate(∅::EmptySet, ::Type{<:HPolygon}, ::Real=Inf; kwargs...) = ∅

overapproximate(∅::EmptySet, ::Real) = ∅

function overapproximate(∅::EmptySet, ::Type{<:Zonotope},
                         ::Union{AbstractDirections,Type{<:AbstractDirections}};
                         kwargs...)
    return ∅
end

for T in (:HPolytope, :HPolyhedron)
    @eval begin
        function overapproximate(∅::EmptySet, ::Type{<:($T)},
                                 dirs::AbstractDirections; prune::Bool=true)
            return ∅
        end
    end
end

overapproximate(∅::EmptySet, dirs::Type{<:AbstractDirections}; kwargs...) = ∅

overapproximate(∅::EmptySet, dirs::AbstractDirections; prune::Bool=true) = ∅

for (T1, T2) in ((:AbstractPolyhedron, :AbstractPolyhedron),
                 (:AbstractPolytope, :AbstractPolytope),
                 (:AbstractPolytope, :AbstractPolyhedron))
    @eval begin
        function overapproximate(cap::Intersection{N,<:($T1),<:($T2)}, dirs::AbstractDirections;
                                 kwargs...) where {N}
            return overapproximate_cap_helper(first(cap), second(cap), dirs; kwargs...)
        end
    end
end

function overapproximate(cap::Intersection{N,<:AbstractPolyhedron,<:AbstractPolytope},
                         dirs::AbstractDirections; kwargs...) where {N}
    return overapproximate_cap_helper(second(cap), first(cap), dirs; kwargs...)
end
