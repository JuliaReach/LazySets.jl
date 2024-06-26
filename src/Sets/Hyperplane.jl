export Hyperplane

"""
    Hyperplane{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a hyperplane of the form ``a⋅x = b``.

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The plane ``y = 0``:

```jldoctest
julia> Hyperplane([0, 1.], 0.)
Hyperplane{Float64, Vector{Float64}}([0.0, 1.0], 0.0)
```
"""
struct Hyperplane{N,VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    a::VN
    b::N

    function Hyperplane(a::VN, b::N) where {N,VN<:AbstractVector{N}}
        @assert !iszero(a) "a hyperplane needs a non-zero normal vector"
        return new{N,VN}(a, b)
    end
end

isoperationtype(::Type{<:Hyperplane}) = false

"""
    normalize(H::Hyperplane{N}, p::Real=N(2)) where {N}

Normalize a hyperplane.

### Input

- `H` -- hyperplane
- `p` -- (optional, default: `2`) norm

### Output

A new hyperplane whose normal direction ``a`` is normalized, i.e., such that
``‖a‖_p = 1`` holds.
"""
function normalize(H::Hyperplane{N}, p::Real=N(2)) where {N}
    a, b = _normalize_halfspace(H, p)
    return Hyperplane(a, b)
end

"""
    constraints_list(H::Hyperplane)

Return the list of constraints of a hyperplane.

### Input

- `H` -- hyperplane

### Output

A list containing two half-spaces.
"""
function constraints_list(H::Hyperplane)
    return _constraints_list_hyperplane(H.a, H.b)
end

# internal helper function
function _constraints_list_hyperplane(a::AbstractVector, b)
    return [HalfSpace(a, b), HalfSpace(-a, -b)]
end

"""
    dim(H::Hyperplane)

Return the dimension of a hyperplane.

### Input

- `H` -- hyperplane

### Output

The ambient dimension of the hyperplane.
"""
function dim(H::Hyperplane)
    return length(H.a)
end

"""
    ρ(d::AbstractVector, H::Hyperplane)

Evaluate the support function of a hyperplane in a given direction.

### Input

- `d` -- direction
- `H` -- hyperplane

### Output

The support function of the hyperplane.
If the set is unbounded in the given direction, the result is `Inf`.
"""
function ρ(d::AbstractVector, H::Hyperplane)
    v, unbounded = _σ_hyperplane_halfspace(d, H.a, H.b; error_unbounded=false,
                                           halfspace=false)
    if unbounded
        N = promote_type(eltype(d), eltype(H))
        return N(Inf)
    end
    return dot(d, v)
end

"""
    σ(d::AbstractVector, H::Hyperplane)

Return a support vector of a hyperplane.

### Input

- `d` -- direction
- `H` -- hyperplane

### Output

A support vector in the given direction, which is only defined in the following
two cases:
1. The direction has norm zero.
2. The direction is the hyperplane's normal direction or its opposite direction.
In all cases, any point on the hyperplane is a solution.
Otherwise this function throws an error.
"""
function σ(d::AbstractVector, H::Hyperplane)
    v, unbounded = _σ_hyperplane_halfspace(d, H.a, H.b; error_unbounded=true,
                                           halfspace=false)
    return v
end

"""
    isbounded(H::Hyperplane)

Check whether a hyperplane is bounded.

### Input

- `H` -- hyperplane

### Output

`true` iff `H` is one-dimensional.
"""
function isbounded(H::Hyperplane)
    return dim(H) == 1
end

"""
    isuniversal(H::Hyperplane, [witness]::Bool=false)

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
function isuniversal(H::Hyperplane, witness::Bool=false)
    if witness
        v = _non_element_halfspace(H.a, H.b)
        return (false, v)
    else
        return false
    end
end

"""
    an_element(H::Hyperplane)

Return some element of a hyperplane.

### Input

- `H` -- hyperplane

### Output

An element on the hyperplane.
"""
function an_element(H::Hyperplane)
    return _an_element_helper_hyperplane(H.a, H.b)
end

"""
    ∈(x::AbstractVector, H::Hyperplane)

Check whether a given point is contained in a hyperplane.

### Input

- `x` -- point/vector
- `H` -- hyperplane

### Output

`true` iff ``x ∈ H``.

### Algorithm

We just check whether ``x`` satisfies ``a⋅x = b``.
"""
function ∈(x::AbstractVector, H::Hyperplane)
    return _isapprox(dot(H.a, x), H.b)
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
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    a = randn(rng, N, dim)
    while iszero(a)
        a = randn(rng, N, dim)
    end
    b = randn(rng, N)
    return Hyperplane(a, b)
end

"""
    isempty(H::Hyperplane)

Check whether a hyperplane is empty.

### Input

- `H` -- hyperplane

### Output

`false`.
"""
function isempty(H::Hyperplane)
    return false
end

"""
    constrained_dimensions(H::Hyperplane)

Return the dimensions in which a hyperplane is constrained.

### Input

- `H` -- hyperplane

### Output

A vector of ascending indices `i` such that the hyperplane is constrained in
dimension `i`.

### Examples

A 2D hyperplane with constraint ``x_1 = 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(H::Hyperplane)
    return nonzero_indices(H.a)
end

"""
```
    _σ_hyperplane_halfspace(d::AbstractVector, a, b;
                            [error_unbounded]::Bool=true,
                            [halfspace]::Bool=false)
```

Return a support vector of a hyperplane ``a⋅x = b`` in direction `d`.

### Input

- `d`         -- direction
- `a`         -- normal direction
- `b`         -- constraint
- `error_unbounded` -- (optional, default: `true`) `true` if an error should be
                 thrown whenever the set is unbounded in the given direction
- `halfspace` -- (optional, default: `false`) `true` if the support vector
                 should be computed for a half-space

### Output

A pair `(v, f)` where `v` is a vector and `f` is a Boolean flag.

The flag `f` is `false` in one of the following cases:
1. The direction has norm zero.
2. The direction is (a multiple of) the hyperplane's normal direction.
3. The direction is (a multiple of) the opposite of the hyperplane's normal
direction and `halfspace` is `false`.
In all these cases, `v` is any point on the hyperplane.

Otherwise, the flag `f` is `true`, the set is unbounded in the given direction,
and `v` is any vector.

If `error_unbounded` is `true` and the set is unbounded in the given direction,
this function throws an error instead of returning.

### Notes

For correctness, consider the
[weak duality of LPs](https://en.wikipedia.org/wiki/Linear_programming#Duality):
If the primal is unbounded, then the dual is infeasible.
Since there is only a single constraint, the feasible set of the dual problem is
``a ⋅ y == d``, ``y ≥ 0`` (with objective function ``b ⋅ y``).
It is easy to see that this problem is infeasible whenever ``a`` is not parallel
to ``d``.
"""
@inline function _σ_hyperplane_halfspace(d::AbstractVector, a, b;
                                         error_unbounded::Bool=true,
                                         halfspace::Bool=false)
    @assert length(d) == length(a) "cannot compute the support vector of a " *
                                   "$(length(a))-dimensional " *
                                   (halfspace ? "halfspace" : "hyperplane") *
                                   " along a vector of length $(length(d))"

    first_nonzero_entry_a = -1
    unbounded = false
    if iszero(d)
        # zero vector
        return (_an_element_helper_hyperplane(a, b), false)
    else
        # not the zero vector, check if it is a normal vector
        N = promote_type(eltype(d), eltype(a))
        factor = zero(N)
        for i in eachindex(a)
            if a[i] == 0
                if d[i] != 0
                    unbounded = true
                    break
                end
            else
                if d[i] == 0
                    unbounded = true
                    break
                elseif first_nonzero_entry_a == -1
                    factor = a[i] / d[i]
                    first_nonzero_entry_a = i
                    if halfspace && factor < 0
                        unbounded = true
                        break
                    end
                elseif d[i] * factor != a[i]
                    unbounded = true
                    break
                end
            end
        end
        if !unbounded
            return (_an_element_helper_hyperplane(a, b, first_nonzero_entry_a), false)
        end
        if error_unbounded
            error("the support vector for the " *
                  (halfspace ? "halfspace" : "hyperplane") * " with normal " *
                  "direction $a is not defined along a direction $d")
        end
        # the first return value does not have a meaning here
        return (d, true)
    end
end

"""
    _an_element_helper_hyperplane(a::AbstractVector{N}, b,
                                  [nonzero_entry_a]::Int) where {N}

Helper function that computes an element on a hyperplane ``a⋅x = b``.

### Input

- `a`               -- normal direction
- `b`               -- constraint
- `nonzero_entry_a` -- (optional, default: computes the first index) index `i`
                       such that `a[i]` is different from 0

### Output

An element on a hyperplane.

### Algorithm

We compute the point on the hyperplane as follows:
- We already found a nonzero entry of ``a`` in dimension, say, ``i``.
- We set ``x[i] = b / a[i]``.
- We set ``x[j] = 0`` for all ``j ≠ i``.
"""
@inline function _an_element_helper_hyperplane(a::AbstractVector{N}, b,
                                               nonzero_entry_a::Int=findfirst(!iszero, a)) where {N}
    x = zeros(N, length(a))
    x[nonzero_entry_a] = b / a[nonzero_entry_a]
    return x
end

function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::Hyperplane{N},
                                 algo::AbstractLinearMapAlgorithm) where {N}
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
    translate(H::Hyperplane, v::AbstractVector; share::Bool=false)

Translate (i.e., shift) a hyperplane by a given vector.

### Input

- `H`     -- hyperplane
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated hyperplane.

### Notes

The normal vector of the hyperplane (vector ``a`` in ``a⋅x = b``) is shared with
the original hyperplane if `share == true`.

### Algorithm

A hyperplane ``a⋅x = b`` is transformed to the hyperplane ``a⋅x = b + a⋅v``.
In other words, we add the dot product ``a⋅v`` to ``b``.
"""
function translate(H::Hyperplane, v::AbstractVector; share::Bool=false)
    @assert length(v) == dim(H) "cannot translate a $(dim(H))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    a = share ? H.a : copy(H.a)
    b = H.b + dot(H.a, v)
    return Hyperplane(a, b)
end

function project(H::Hyperplane{N}, block::AbstractVector{Int}; kwargs...) where {N}
    if constrained_dimensions(H) ⊆ block
        return Hyperplane(H.a[block], H.b)
    else
        return Universe{N}(length(block))
    end
end

"""
    project(x::AbstractVector, H::Hyperplane)

Project a point onto a hyperplane.

### Input

- `x` -- point
- `H` -- hyperplane

### Output

The projection of `x` onto `H`.

### Algorithm

The projection of ``x`` onto the hyperplane of the form ``a⋅x = b`` is

```math
    x - \\dfrac{a (a⋅x - b)}{‖a‖²}
```
"""
function project(x::AbstractVector, H::Hyperplane)
    return x - H.a * (dot(H.a, x) - H.b) / norm(H.a, 2)^2
end

function is_hyperplanar(::Hyperplane)
    return true
end

# ============================================
# Functionality that requires Symbolics
# ============================================
function load_symbolics_hyperplane()
    return quote

        # returns `(true, sexpr)` if expr represents a hyperplane,
        # where sexpr is the simplified expression sexpr := LHS - RHS == 0
        # otherwise returns `(false, expr)`
        function _is_hyperplane(expr::Symbolic)
            got_hyperplane = operation(expr) == ==
            if got_hyperplane
                # simplify to the form a*x + b == 0
                a, b = arguments(expr)
                sexpr = simplify(a - b)
            end
            return got_hyperplane ? (true, sexpr) : (false, expr)
        end

        """
            Hyperplane(expr::Num, vars=_get_variables(expr); [N]::Type{<:Real}=Float64)

        Return the hyperplane given by a symbolic expression.

        ### Input

        - `expr` -- symbolic expression that describes a hyperplane
        - `vars` -- (optional, default: `_get_variables(expr)`), if a vector of
                    variables is given, use those as the ambient variables with respect
                    to which derivations take place; otherwise, use only the variables
                    that appear in the given expression (but be careful because the
                    order may be incorrect; it is advised to always specify `vars`
                    explicitly)
        - `N`    -- (optional, default: `Float64`) the numeric type of the hyperplane

        ### Output

        A `Hyperplane`.

        ### Examples

        ```jldoctest
        julia> using Symbolics

        julia> vars = @variables x y
        2-element Vector{Num}:
         x
         y

        julia> Hyperplane(x - y == 2)
        Hyperplane{Float64, Vector{Float64}}([-1.0, 1.0], 2.0)

        julia> Hyperplane(x == y)
        Hyperplane{Float64, Vector{Float64}}([1.0, -1.0], -0.0)

        julia> vars = @variables x[1:4]
        1-element Vector{Symbolics.Arr{Num, 1}}:
         x[1:4]

        julia> Hyperplane(x[1] == x[2], x)
        Hyperplane{Float64, Vector{Float64}}([1.0, -1.0, 0.0, 0.0], -0.0)
        ```

        ### Algorithm

        It is assumed that the expression is of the form
        `α*x1 + ⋯ + α*xn + γ == β*x1 + ⋯ + β*xn + δ`.
        This expression is transformed, by rearrangement and substitution, into the
        canonical form `a1 * x1 + ⋯ + an * xn == b`. To identify the coefficients, we
        take derivatives with respect to the ambient variables `vars`. Therefore, the
        order in which the variables appear in `vars` affects the final result. Finally,
        the returned set is the hyperplane with normal vector `[a1, …, an]` and
        displacement `b`.
        """
        function Hyperplane(expr::Num, vars::AbstractVector{Num}=_get_variables(expr);
                            N::Type{<:Real}=Float64)
            valid, sexpr = _is_hyperplane(Symbolics.value(expr))
            if !valid
                throw(ArgumentError("expected an expression of the form `ax == b`, got $expr"))
            end

            # compute the linear coefficients by taking first order derivatives
            coeffs = [N(α.val) for α in gradient(sexpr, collect(vars))]

            # get the constant term by expression substitution
            zeroed_vars = Dict(v => zero(N) for v in vars)
            β = -N(Symbolics.substitute(sexpr, zeroed_vars))

            return Hyperplane(coeffs, β)
        end

        Hyperplane(expr::Num, vars; N::Type{<:Real}=Float64) = Hyperplane(expr, _vec(vars); N=N)
    end
end  # quote / load_symbolics_hyperplane()

# =====================================
# Functionality that requires SymEngine
# =====================================

function load_symengine_hyperplane()
    return quote
        """
            _is_hyperplane(expr::Expr)

        Determine whether the given expression corresponds to a hyperplane.

        ### Input

        - `expr` -- a symbolic expression

        ### Output

        `true` if `expr` corresponds to a half-space or `false` otherwise.

        ### Examples

        ```jldoctest
        julia> using LazySets: _is_hyperplane

        julia> _is_hyperplane(:(x1 = 0))
        true

        julia> _is_hyperplane(:(2*x1 = 4))
        true

        julia> _is_hyperplane(:(6.1 = 5.3*f - 0.1*g))
        true

        julia> _is_hyperplane(:(2*x1^2 = 4))
        false

        julia> _is_hyperplane(:(x1^2 = 4*x2 - x3))
        false

        julia> _is_hyperplane(:(x1 = 4*x2 - x3))
        true
        ```
        """
        function _is_hyperplane(expr::Expr)::Bool

            # check that there are three arguments
            # these are the comparison symbol, the left hand side and the right hand side
            if (length(expr.args) != 2) || !(expr.head == :(=))
                return false
            end

            # convert to symengine expressions
            lhs = convert(Basic, expr.args[1])

            if :args in fieldnames(typeof(expr.args[2]))
                # treats the 4 in :(2*x1 = 4)
                rhs = convert(Basic, expr.args[2].args[2])
            else
                rhs = convert(Basic, expr.args[2])
            end

            # check if the expression defines a hyperplane
            return _is_linearcombination(lhs) && _is_linearcombination(rhs)
        end

        """
            convert(::Type{Hyperplane{N}}, expr::Expr; vars=nothing) where {N}

        Return a `LazySet.Hyperplane` given a symbolic expression that represents a hyperplane.

        ### Input

        - `expr` -- a symbolic expression
        - `vars` -- (optional, default: `nothing`): set of variables with respect to which
                    the gradient is taken; if nothing, it takes the free symbols in the given expression

        ### Output

        A `Hyperplane`, in the form `ax = b`.

        ### Examples

        ```jldoctest convert_hyperplane
        julia> convert(Hyperplane, :(x1 = -0.03))
        Hyperplane{Float64, Vector{Float64}}([1.0], -0.03)

        julia> convert(Hyperplane, :(x1 + 0.03 = 0))
        Hyperplane{Float64, Vector{Float64}}([1.0], -0.03)

        julia> convert(Hyperplane, :(x1 + x2 = 2*x4 + 6))
        Hyperplane{Float64, Vector{Float64}}([1.0, 1.0, -2.0], 6.0)
        ```

        You can also specify the set of "ambient" variables in the hyperplane, even if not
        all of them appear:

        ```jldoctest convert_hyperplane
        julia> using SymEngine: Basic

        julia> convert(Hyperplane, :(x1 + x2 = 2*x4 + 6), vars=Basic[:x1, :x2, :x3, :x4])
        Hyperplane{Float64, Vector{Float64}}([1.0, 1.0, 0.0, -2.0], 6.0)
        ```
        """
        function convert(::Type{Hyperplane{N}}, expr::Expr; vars::Vector{Basic}=Basic[]) where {N}
            @assert _is_hyperplane(expr) "the expression :(expr) does not correspond to a Hyperplane"

            # get sides of the inequality
            lhs = convert(Basic, expr.args[1])

            # treats the 4 in :(2*x1 = 4)
            rhs = :args in fieldnames(typeof(expr.args[2])) ? convert(Basic, expr.args[2].args[2]) :
                convert(Basic, expr.args[2])

            # a1 x1 + ... + an xn + K = 0
            eq = lhs - rhs
            if isempty(vars)
                vars = free_symbols(eq)
            end
            K = subs(eq, [vi => zero(N) for vi in vars]...)
            a = convert(Basic, eq - K)

            # convert to numeric types
            K = convert(N, K)
            a = convert(Vector{N}, diff.(a, vars))

            return Hyperplane(a, -K)
        end

        # type-less default Hyperplane conversion
        function convert(::Type{Hyperplane}, expr::Expr; vars::Vector{Basic}=Basic[])
            return convert(Hyperplane{Float64}, expr; vars=vars)
        end

        function free_symbols(expr::Expr, ::Type{<:Hyperplane})
            # get sides of the equality
            lhs = convert(Basic, expr.args[1])

            # treats the 4 in :(2*x1 = 4)
            rhs = :args in fieldnames(typeof(expr.args[2])) ? convert(Basic, expr.args[2].args[2]) :
                  convert(Basic, expr.args[2])
            return free_symbols(lhs - rhs)
        end
    end
end  # quote / load_symengine_hyperplane()

"""
    distance(x::AbstractVector, H::Hyperplane{N}) where {N}

Compute the distance between point `x` and hyperplane `H` with respect to the
Euclidean norm.

### Input

- `x` -- vector
- `H` -- hyperplane

### Output

A scalar representing the distance between point `x` and hyperplane `H`.
"""
@commutative function distance(x::AbstractVector, H::Hyperplane{N}) where {N}
    a, b = _normalize_halfspace(H, N(2))
    return abs(dot(x, a) - b)
end

"""
    reflect(x::AbstractVector, H::Hyperplane)

Reflect (mirror) a vector in a hyperplane.

### Input

- `x` -- point/vector
- `H` -- hyperplane

### Output

The reflection of `x` in `H`.

### Algorithm

The reflection of a point ``x`` in the hyperplane ``a ⋅ x = b`` is

```math
    x − 2 \\frac{x ⋅ a − b}{a ⋅ a} a
```

where ``u · v`` denotes the dot product.
"""
@commutative function reflect(x::AbstractVector, H::Hyperplane)
    return _reflect_point_hyperplane(x, H.a, H.b)
end

function _reflect_point_hyperplane(x, a, b)
    return x - 2 * (dot(x, a) - b) / dot(a, a) * a
end
