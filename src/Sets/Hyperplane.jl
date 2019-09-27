import Base: rand,
             ∈,
             isempty

export Hyperplane,
       an_element

"""
    Hyperplane{N<:Real} <: AbstractPolyhedron{N}

Type that represents a hyperplane of the form ``a⋅x = b``.

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The plane ``y = 0``:

```jldoctest
julia> Hyperplane([0, 1.], 0.)
Hyperplane{Float64}([0.0, 1.0], 0.0)
```
"""
struct Hyperplane{N<:Real} <: AbstractPolyhedron{N}
    a::AbstractVector{N}
    b::N

    function Hyperplane{N}(a::AbstractVector{N}, b::N) where {N<:Real}
        @assert !iszero(a) "a hyperplane needs a non-zero normal vector"
        return new{N}(a, b)
    end
end

isoperationtype(::Type{<:Hyperplane}) = false

# convenience constructor without type parameter
Hyperplane(a::AbstractVector{N}, b::N) where {N<:Real} = Hyperplane{N}(a, b)


# --- polyhedron interface functions ---


"""
    constraints_list(hp::Hyperplane{N}) where {N<:Real}

Return the list of constraints of a hyperplane.

### Input

- `hp` -- hyperplane

### Output

A list containing two half-spaces.
"""
function constraints_list(hp::Hyperplane{N}) where {N<:Real}
    return _constraints_list_hyperplane(hp.a, hp.b)
end


# --- LazySet interface functions ---


"""
    dim(hp::Hyperplane)::Int

Return the dimension of a hyperplane.

### Input

- `hp` -- hyperplane

### Output

The ambient dimension of the hyperplane.
"""
function dim(hp::Hyperplane)::Int
    return length(hp.a)
end

"""
    ρ(d::AbstractVector{N}, hp::Hyperplane{N})::N where {N<:Real}

Evaluate the support function of a hyperplane in a given direction.

### Input

- `d`  -- direction
- `hp` -- hyperplane

### Output

The support function of the hyperplane.
If the set is unbounded in the given direction, the result is `Inf`.
"""
function ρ(d::AbstractVector{N}, hp::Hyperplane{N})::N where {N<:Real}
    v, unbounded = σ_helper(d, hp, error_unbounded=false)
    if unbounded
        return N(Inf)
    end
    return dot(d, v)
end

"""
    σ(d::AbstractVector{N}, hp::Hyperplane{N}) where {N<:Real}

Return the support vector of a hyperplane.

### Input

- `d`  -- direction
- `hp` -- hyperplane

### Output

The support vector in the given direction, which is only defined in the
following two cases:
1. The direction has norm zero.
2. The direction is the hyperplane's normal direction or its opposite direction.
In all cases, the result is any point on the hyperplane.
Otherwise this function throws an error.
"""
function σ(d::AbstractVector{N}, hp::Hyperplane{N}) where {N<:Real}
    v, unbounded = σ_helper(d, hp, error_unbounded=true)
    return v
end

"""
    isbounded(hp::Hyperplane)::Bool

Determine whether a hyperplane is bounded.

### Input

- `hp` -- hyperplane

### Output

`false`.
"""
function isbounded(::Hyperplane)::Bool
    return false
end

"""
    isuniversal(hp::Hyperplane{N}, [witness]::Bool=false
               )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a hyperplane is universal.

### Input

- `P`       -- hyperplane
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ P``

### Algorithm

A witness is produced by adding the normal vector to an element on the
hyperplane.
"""
function isuniversal(hp::Hyperplane{N}, witness::Bool=false
                    )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    if witness
        v = an_element(hp) + hp.a
        return (false, v)
    else
        return false
    end
end

"""
    an_element(hp::Hyperplane{N})::Vector{N} where {N<:Real}

Return some element of a hyperplane.

### Input

- `hp` -- hyperplane

### Output

An element on the hyperplane.
"""
function an_element(hp::Hyperplane{N})::Vector{N} where {N<:Real}
    return an_element_helper(hp)
end

"""
    ∈(x::AbstractVector{N}, hp::Hyperplane{N})::Bool where {N<:Real}

Check whether a given point is contained in a hyperplane.

### Input

- `x` -- point/vector
- `hp` -- hyperplane

### Output

`true` iff ``x ∈ hp``.

### Algorithm

We just check if ``x`` satisfies ``a⋅x = b``.
"""
function ∈(x::AbstractVector{N}, hp::Hyperplane{N})::Bool where {N<:Real}
    return dot(x, hp.a) == hp.b
end

"""
    rand(::Type{Hyperplane}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::Hyperplane{N}

Create a random hyperplane.

### Input

- `Hyperplane` -- type for dispatch
- `N`          -- (optional, default: `Float64`) numeric type
- `dim`        -- (optional, default: 2) dimension
- `rng`        -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`       -- (optional, default: `nothing`) seed for reseeding

### Output

A random hyperplane.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the constraint `a` is nonzero.
"""
function rand(::Type{Hyperplane};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::Hyperplane{N}
    rng = reseed(rng, seed)
    a = randn(rng, N, dim)
    while iszero(a)
        a = randn(rng, N, dim)
    end
    b = randn(rng, N)
    return Hyperplane(a, b)
end

"""
    isempty(hp::Hyperplane)::Bool

Return if a hyperplane is empty or not.

### Input

- `hp` -- hyperplane

### Output

`false`.
"""
function isempty(hp::Hyperplane)::Bool
    return false
end

"""
    constrained_dimensions(hp::Hyperplane{N})::Vector{Int} where {N<:Real}

Return the indices in which a hyperplane is constrained.

### Input

- `hp` -- hyperplane

### Output

A vector of ascending indices `i` such that the hyperplane is constrained in
dimension `i`.

### Examples

A 2D hyperplane with constraint ``x1 = 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(hp::Hyperplane{N})::Vector{Int} where {N<:Real}
    return nonzero_indices(hp.a)
end


# --- Hyperplane functions ---


"""
```
    σ_helper(d::AbstractVector{N},
             hp::Hyperplane{N};
             error_unbounded::Bool=true,
             [halfspace]::Bool=false) where {N<:Real}
```

Return the support vector of a hyperplane.

### Input

- `d`         -- direction
- `hp`        -- hyperplane
- `error_unbounded` -- (optional, default: `true`) `true` if an error should be
                 thrown whenever the set is
                 unbounded in the given direction
- `halfspace` -- (optional, default: `false`) `true` if the support vector
                 should be computed for a half-space

### Output

A pair `(v, b)` where `v` is a vector and `b` is a Boolean flag.

The flag `b` is `false` in one of the following cases:
1. The direction has norm zero.
2. The direction is the hyperplane's normal direction.
3. The direction is the opposite of the hyperplane's normal direction and
`halfspace` is `false`.
In all these cases, `v` is any point on the hyperplane.

Otherwise, the flag `b` is `true`, the set is unbounded in the given direction,
and `v` is any vector.

If `error_unbounded` is `true` and the set is unbounded in the given direction,
this function throws an error instead of returning.

### Notes

For correctness, consider the [weak duality of
LPs](https://en.wikipedia.org/wiki/Linear_programming#Duality):
If the primal is unbounded, then the dual is infeasible.
Since there is only a single constraint, the feasible set of the dual problem is
`hp.a ⋅ y == d`, `y >= 0` (with objective function `hp.b ⋅ y`).
It is easy to see that this problem is infeasible whenever `a` is not parallel
to `d`.
"""
@inline function σ_helper(d::AbstractVector{N},
                          hp::Hyperplane{N};
                          error_unbounded::Bool=true,
                          halfspace::Bool=false) where {N<:Real}
    @assert (length(d) == dim(hp)) "cannot compute the support vector of a " *
        "$(dim(hp))-dimensional " * (halfspace ? "halfspace" : "hyperplane") *
        " along a vector of length $(length(d))"

    first_nonzero_entry_a = -1
    unbounded = false
    if iszero(d)
        # zero vector
        return (an_element(hp), false)
    else
        # not the zero vector, check if it is a normal vector
        factor = zero(N)
        for i in 1:length(hp.a)
            if hp.a[i] == 0
                if d[i] != 0
                    unbounded = true
                    break
                end
            else
                if d[i] == 0
                    unbounded = true
                    break
                elseif first_nonzero_entry_a == -1
                    factor = hp.a[i] / d[i]
                    first_nonzero_entry_a = i
                    if halfspace && factor < 0
                        unbounded = true
                        break
                    end
                elseif d[i] * factor != hp.a[i]
                    unbounded = true
                    break
                end
            end
        end
        if !unbounded
            return (an_element_helper(hp, first_nonzero_entry_a), false)
        end
        if error_unbounded
            error("the support vector for the " *
                (halfspace ? "halfspace" : "hyperplane") * " with normal " *
                "direction $(hp.a) is not defined along a direction $d")
        end
        # the first return value does not have a meaning here
        return (d, true)
    end
end

"""
    an_element_helper(hp::Hyperplane{N},
                      [nonzero_entry_a]::Int)::Vector{N} where {N<:Real}

Helper function that computes an element on a hyperplane's hyperplane.

### Input

- `hp` -- hyperplane
- `nonzero_entry_a` -- (optional, default: computes the first index) index `i`
                       such that `hp.a[i]` is different from 0

### Output

An element on a hyperplane.

### Algorithm

We compute the point on the hyperplane as follows:
- We already found a nonzero entry of ``a`` in dimension, say, ``i``.
- We set ``x[i] = b / a[i]``.
- We set ``x[j] = 0`` for all ``j ≠ i``.
"""
@inline function an_element_helper(hp::Hyperplane{N},
                                   nonzero_entry_a::Int=findnext(x -> x!=zero(N), hp.a, 1)
                                  )::Vector{N} where {N<:Real}
    @assert nonzero_entry_a in 1:length(hp.a) "invalid index " *
        "$nonzero_entry_a for hyperplane"
    x = zeros(N, dim(hp))
    x[nonzero_entry_a] = hp.b / hp.a[nonzero_entry_a]
    return x
end

# internal helper function
function _constraints_list_hyperplane(a::AbstractVector{N}, b::N
                                     ) where {N<:Real}
    return [HalfSpace(a, b), HalfSpace(-a, -b)]
end

function _linear_map_hrep(M::AbstractMatrix{N}, P::Hyperplane{N}, use_inv::Bool
                         ) where {N<:Real}
    constraint = _linear_map_hrep_helper(M, P, use_inv)[1]
    return Hyperplane(constraint.a, constraint.b)
end

"""
    translate(hp::Hyperplane{N}, v::AbstractVector{N}; share::Bool=false
             ) where {N<:Real}

Translate (i.e., shift) a hyperplane by a given vector.

### Input

- `hp`    -- hyperplane
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated hyperplane.

### Notes

The normal vectors of the hyperplane (vector `a` in `a⋅x = b`) is shared with
the original hyperplane if `share == true`.

### Algorithm

A hyperplane ``a⋅x = b`` is transformed to the hyperplane ``a⋅x = b + a⋅v``.
In other words, we add the dot product ``a⋅v`` to ``b``.
"""
function translate(hp::Hyperplane{N}, v::AbstractVector{N}; share::Bool=false
                  ) where {N<:Real}
    @assert length(v) == dim(hp) "cannot translate a $(dim(hp))-dimensional " *
                                 "set by a $(length(v))-dimensional vector"
    a = share ? hp.a : copy(hp.a)
    b = hp.b + dot(hp.a, v)
    return Hyperplane(a, b)
end
