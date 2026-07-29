using Symbolics: Num, arguments, gradient, operation, simplify, substitute, value
import LazySets.HyperplaneModule: Hyperplane

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

# returns `(true, sexpr)` if expr represents a hyperplane,
# where sexpr is the simplified expression sexpr := LHS - RHS == 0
# otherwise returns `(false, expr)`
function _ishyperplanar(expr::BasicSymbolic2)
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
julia> using LazySets, Symbolics

julia> vars = @variables x y
2-element Vector{Num}:
 x
 y

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
    valid, sexpr = _ishyperplanar(value(expr))
    if !valid
        throw(ArgumentError("expected an expression of the form `ax == b`, got $expr"))
    end

    # compute the linear coefficients by taking first order derivatives
    coeffs = [N(value(α)) for α in gradient(sexpr, collect(vars))]

    # get the constant term by expression substitution
    zeroed_vars = Dict(v => zero(N) for v in vars)
    β = -N(value(substitute(sexpr, zeroed_vars)))

    return Hyperplane(coeffs, β)
end

function Hyperplane(expr::Num, vars; N::Type{<:Real}=Float64)
    return Hyperplane(expr, _vec(vars); N=N)
end
