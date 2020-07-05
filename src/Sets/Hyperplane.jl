import Base: rand,
             ∈,
             isempty

export Hyperplane,
       an_element

"""
    Hyperplane{N<:Real, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a hyperplane of the form ``a⋅x = b``.

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The plane ``y = 0``:

```jldoctest
julia> Hyperplane([0, 1.], 0.)
Hyperplane{Float64,Array{Float64,1}}([0.0, 1.0], 0.0)
```
"""
struct Hyperplane{N<:Real, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    a::VN
    b::N

    function Hyperplane(a::VN, b::N) where {N<:Real, VN<:AbstractVector{N}}
        @assert !iszero(a) "a hyperplane needs a non-zero normal vector"
        return new{N, VN}(a, b)
    end
end

isoperationtype(::Type{<:Hyperplane}) = false
isconvextype(::Type{<:Hyperplane}) = true


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
    dim(hp::Hyperplane)

Return the dimension of a hyperplane.

### Input

- `hp` -- hyperplane

### Output

The ambient dimension of the hyperplane.
"""
function dim(hp::Hyperplane)
    return length(hp.a)
end

"""
    ρ(d::AbstractVector{N}, hp::Hyperplane{N}) where {N<:Real}

Evaluate the support function of a hyperplane in a given direction.

### Input

- `d`  -- direction
- `hp` -- hyperplane

### Output

The support function of the hyperplane.
If the set is unbounded in the given direction, the result is `Inf`.
"""
function ρ(d::AbstractVector{N}, hp::Hyperplane{N}) where {N<:Real}
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
    isbounded(hp::Hyperplane)

Determine whether a hyperplane is bounded.

### Input

- `hp` -- hyperplane

### Output

`false`.
"""
function isbounded(::Hyperplane)
    return false
end

"""
    isuniversal(hp::Hyperplane{N}, [witness]::Bool=false) where {N<:Real}

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
function isuniversal(hp::Hyperplane{N}, witness::Bool=false) where {N<:Real}
    if witness
        v = an_element(hp) + hp.a
        return (false, v)
    else
        return false
    end
end

"""
    an_element(hp::Hyperplane{N}) where {N<:Real}

Return some element of a hyperplane.

### Input

- `hp` -- hyperplane

### Output

An element on the hyperplane.
"""
function an_element(hp::Hyperplane{N}) where {N<:Real}
    return an_element_helper(hp)
end

"""
    ∈(x::AbstractVector{N}, hp::Hyperplane{N}) where {N<:Real}

Check whether a given point is contained in a hyperplane.

### Input

- `x` -- point/vector
- `hp` -- hyperplane

### Output

`true` iff ``x ∈ hp``.

### Algorithm

We just check if ``x`` satisfies ``a⋅x = b``.
"""
function ∈(x::AbstractVector{N}, hp::Hyperplane{N}) where {N<:Real}
    return dot(x, hp.a) == hp.b
end

"""
    rand(::Type{Hyperplane}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

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
              seed::Union{Int, Nothing}=nothing)
    rng = reseed(rng, seed)
    a = randn(rng, N, dim)
    while iszero(a)
        a = randn(rng, N, dim)
    end
    b = randn(rng, N)
    return Hyperplane(a, b)
end

"""
    isempty(hp::Hyperplane)

Return if a hyperplane is empty or not.

### Input

- `hp` -- hyperplane

### Output

`false`.
"""
function isempty(hp::Hyperplane)
    return false
end

"""
    constrained_dimensions(hp::Hyperplane{N}) where {N<:Real}

Return the indices in which a hyperplane is constrained.

### Input

- `hp` -- hyperplane

### Output

A vector of ascending indices `i` such that the hyperplane is constrained in
dimension `i`.

### Examples

A 2D hyperplane with constraint ``x1 = 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(hp::Hyperplane{N}) where {N<:Real}
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
                      [nonzero_entry_a]::Int) where {N<:Real}

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
                                  ) where {N<:Real}
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

function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::Hyperplane{N},
                                 algo::AbstractLinearMapAlgorithm) where {N<:Real}
    constraints = _linear_map_hrep(M, P, algo)
    if length(constraints) == 2
        # assuming these constraints define a hyperplane
        c = first(constraints)
        return Hyperplane(c.a, c.b)
    elseif isempty(constraints)
        return Universe{N}(size(M, 1))
    else
        error("unexpected number of $(length(constraints)) constraints")
    end
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

# ============================================
# Functionality that requires ModelingToolkit
# ============================================
function load_modeling_toolkit_hyperplane()
return quote

"""
    Hyperplane(expr::Operation, vars=get_variables(expr); N::Type{<:Real}=Float64)

Return the hyperplane given by a symbolic expression.

### Input

- `expr` -- symbolic expression that describes a hyperplane
- `vars` -- (optional, default: `get_variables(expr)`), if an array of variables is given,
            use those as the ambient variables in the set with respect to which derivations
            take place; otherwise, use only the variables which appear in the given
            expression (but be careful because the order may change; in doubt, pass `vars` explicitly)
- `N`    -- (optional, default: `Float64`) the numeric type of the returned hyperplane

### Output

A `Hyperplane`.

### Examples

```julia
julia> using ModelingToolkit

julia> vars = @variables x y
(x, y)

julia> Hyperplane(x - y == 2)
Hyperplane{Float64,Array{Float64,1}}([1.0, -1.0], 2.0)

julia> Hyperplane(x == y)
Hyperplane{Float64,Array{Float64,1}}([1.0, -1.0], -0.0)

julia> vars = @variables x[1:4]
(Operation[x₁, x₂, x₃, x₄],)

julia> Hyperplane(x[1] == x[2], x)
Hyperplane{Float64,Array{Float64,1}}([1.0, -1.0, 0.0, 0.0], -0.0)
```

### Algorithm

It is assumed that the expression is of the form
`EXPR0: α*x1 + ⋯ + α*xn + γ == β*x1 + ⋯ + β*xn + δ`.
This expression is transformed, by rearrangement and substitution, into the
canonical form `EXPR1 : a1 * x1 + ⋯ + an * xn == b`. The method used to identify
the coefficients is to take derivatives with respect to the ambient variables `vars`.
Therefore, the order in which the variables appear in `vars` affects the final result.
Finally, the returned set is the hyperplane with normal vector `[a1, …, an]` and
displacement `b`.
"""
function Hyperplane(expr::Operation, vars=get_variables(expr); N::Type{<:Real}=Float64)
    (expr.op == ==) || throw(ArgumentError("expected an expression of the form `ax == b`, got $expr"))

    # simplify to the form a*x + β == 0
    a, b = expr.args
    sexpr = simplify(a - b)

    # compute the linear coefficients by taking first order derivatives
    coeffs = [N(α.value) for α in gradient(sexpr, collect(vars))]

    # get the constant term by expression substitution
    dvars = Dict(to_symbolic(vi) => zero(N) for vi in vars)
    β = -N(ModelingToolkit.SymbolicUtils.substitute(to_symbolic(sexpr), dvars, fold=true))

    return Hyperplane(coeffs, β)
end

end end  # quote / load_modeling_toolkit_hyperplane()
