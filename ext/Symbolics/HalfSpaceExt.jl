using Symbolics: Num, arguments, gradient, operation, simplify, substitute, value
import LazySets.HalfSpaceModule: HalfSpace

if isdefined(SymbolicUtils, :Symbolic)
    using SymbolicUtils: Symbolic
    const BasicSymbolic2 = Symbolic
elseif isdefined(SymbolicUtils, :BasicSymbolic)
    using SymbolicUtils: BasicSymbolic
    const BasicSymbolic2 = BasicSymbolic
else
    using Symbolics: BasicSymbolic
    const BasicSymbolic2 = BasicSymbolic
end

# returns `(true, sexpr)` if expr represents a half-space,
# where sexpr is the simplified expression sexpr := LHS - RHS <= 0
# otherwise, returns `(false, expr)`
function _ishalfspace(expr::BasicSymbolic2)
    got_halfspace = true

    # find sense and normalize
    op = operation(expr)
    args = arguments(expr)
    if op in (<=, <)
        a, b = args
        sexpr = simplify(a - b)

    elseif op in (>=, >)
        a, b = args
        sexpr = simplify(b - a)

    elseif (op == |) && (operation(args[1]) == <)
        a, b = arguments(args[1])
        sexpr = simplify(a - b)

    elseif (op == |) && (operation(args[2]) == <)
        a, b = arguments(args[2])
        sexpr = simplify(a - b)

    elseif (op == |) && (operation(args[1]) == >)
        a, b = arguments(args[1])
        sexpr = simplify(b - a)

    elseif (op == |) && (operation(args[2]) == >)
        a, b = arguments(args[2])
        sexpr = simplify(b - a)

    else
        got_halfspace = false
    end

    return got_halfspace ? (true, sexpr) : (false, expr)
end

"""
    HalfSpace(expr::Num, vars=_get_variables(expr); [N]::Type{<:Real}=Float64)

Return the half-space given by a symbolic expression.

### Input

- `expr` -- symbolic expression that describes a half-space
- `vars` -- (optional, default: `_get_variables(expr)`) if an array of variables
            is given, use those as the ambient variables in the set with respect
            to which derivations take place; otherwise, use only the variables
            that appear in the given expression (but be careful because the
            order may be incorrect; it is advised to always pass `vars`
            explicitly; see the examples below for details)
- `N`    -- (optional, default: `Float64`) the numeric type of the half-space

### Output

A `HalfSpace`.

### Examples

```jldoctest halfspace_symbolics
julia> using LazySets, Symbolics

julia> vars = @variables x y
2-element Vector{Num}:
 x
 y

julia> HalfSpace(x - y <= 2, vars)
HalfSpace{Float64, Vector{Float64}}([1.0, -1.0], 2.0)

julia> HalfSpace(x >= y, vars)
HalfSpace{Float64, Vector{Float64}}([-1.0, 1.0], -0.0)

julia> vars = @variables x[1:4]
1-element Vector{Symbolics.Arr{Num, 1}}:
 x[1:4]

julia> HalfSpace(x[1] >= x[2], x)
HalfSpace{Float64, Vector{Float64}}([-1.0, 1.0, 0.0, 0.0], -0.0)
```

Be careful with using the default `vars` value because it may introduce a wrong
order.

```@example halfspace_symbolics  # doctest deactivated due to downgrade of Symbolics
julia> HalfSpace(2x ≥ 5y - 1) # correct
HalfSpace{Float64, Vector{Float64}}([-2.0, 5.0], 1.0)

julia> HalfSpace(2x ≥ 5y - 1, vars) # correct
HalfSpace{Float64, Vector{Float64}}([-2.0, 5.0], 1.0)

julia> HalfSpace(y - x ≥ 1) # incorrect
HalfSpace{Float64, Vector{Float64}}([-1.0, 1.0], -1.0)

julia> HalfSpace(y - x ≥ 1, vars) # correct
HalfSpace{Float64, Vector{Float64}}([1.0, -1.0], -1.0)
julia> nothing  # hide
```

### Algorithm

It is assumed that the expression is of the form
`α*x1 + ⋯ + α*xn + γ CMP β*x1 + ⋯ + β*xn + δ`,
where `CMP` is one among `<`, `<=`, `≤`, `>`, `>=` or `≥`.
This expression is transformed, by rearrangement and substitution, into the
canonical form `a1 * x1 + ⋯ + an * xn ≤ b`. The method used to identify the
coefficients is to take derivatives with respect to the ambient variables `vars`.
Therefore, the order in which the variables appear in `vars` affects the final
result. Note in particular that strict inequalities are relaxed as being
smaller-or-equal. Finally, the returned set is the half-space with normal vector
`[a1, …, an]` and displacement `b`.
"""
function HalfSpace(expr::Num, vars::AbstractVector{Num}=_get_variables(expr);
                   N::Type{<:Real}=Float64)
    valid, sexpr = _ishalfspace(value(expr))
    if !valid
        throw(ArgumentError("expected an expression describing a half-space, got $expr"))
    end

    # compute the linear coefficients by taking first-order derivatives
    coeffs = [N(value(α)) for α in gradient(sexpr, collect(vars))]

    # get the constant term by expression substitution
    zeroed_vars = Dict(v => zero(N) for v in vars)
    β = -N(value(substitute(sexpr, zeroed_vars)))

    return HalfSpace(coeffs, β)
end

function HalfSpace(expr::Num, vars; N::Type{<:Real}=Float64)
    return HalfSpace(expr, _vec(vars); N=N)
end
