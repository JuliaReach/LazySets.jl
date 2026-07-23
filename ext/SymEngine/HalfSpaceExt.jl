using LazySets.HalfSpaceModule: HalfSpace
using SymEngine: Basic
import Base: convert
import SymEngine: free_symbols

function _parse_halfspace(expr::Expr)
    lhs = convert(Basic, expr.args[2])
    rhs = convert(Basic, expr.args[3])
    cmp = expr.args[1]
    return (lhs - rhs, cmp)
end

function free_symbols(expr::Expr, ::Type{<:HalfSpace})
    linexpr, _ = _parse_halfspace(expr)
    return SymEngine.free_symbols(linexpr)
end

"""
    _ishalfspace(expr::Expr)

Determine whether the given expression corresponds to a half-space.

### Input

- `expr` -- a symbolic expression

### Output

`true` if `expr` corresponds to a half-space or `false` otherwise.

### Examples

```jldoctest
julia> using LazySets

julia> import SymEngine

julia> LazySetsSymEngineExt = Base.get_extension(LazySets, :LazySetsSymEngineExt);

julia> using .LazySetsSymEngineExt: _ishalfspace

julia> all(_ishalfspace.([:(x1 <= 0), :(x1 < 0), :(x1 > 0), :(x1 >= 0)]))
true

julia> _ishalfspace(:(x1 = 0))
false

julia> _ishalfspace(:(2*x1 <= 4))
true

julia> _ishalfspace(:(6.1 <= 5.3*f - 0.1*g))
true

julia> _ishalfspace(:(2*x1^2 <= 4))
false

julia> _ishalfspace(:(x1^2 > 4*x2 - x3))
false

julia> _ishalfspace(:(x1 > 4*x2 - x3))
true
```
"""
function _ishalfspace(expr::Expr)::Bool
    # check that there are three arguments:
    # the comparison symbol, the left-hand side and the right-hand side
    if (length(expr.args) != 3) || !(expr.head == :call)
        return false
    end

    # convert to SymEngine expression
    linexpr, cmp = _parse_halfspace(expr)

    # check that this is an inequality
    if cmp ∉ [:(<=), :(<), :(>=), :(>)]
        return false
    end

    # check if the expression defines a half-space
    return _is_linear_combination(linexpr)
end

"""
    convert(::Type{HalfSpace{N}}, expr::Expr; vars::Vector{Basic}=Basic[]) where {N}

Return a `LazySet.HalfSpace` given a symbolic expression that represents a half-space.

### Input

- `expr` -- a symbolic expression
- `vars` -- (optional, default: `nothing`): set of variables with respect to which the
            gradient is taken; if empty, we take the free symbols in the given expression

### Output

A `HalfSpace`, in the form `ax <= b`.

### Examples

```jldoctest convert_halfspace
julia> using LazySets

julia> import SymEngine

julia> convert(HalfSpace, :(x1 <= -0.03))
HalfSpace{Float64, Vector{Float64}}([1.0], -0.03)

julia> convert(HalfSpace, :(x1 < -0.03))
HalfSpace{Float64, Vector{Float64}}([1.0], -0.03)

julia> convert(HalfSpace, :(x1 > -0.03))
HalfSpace{Float64, Vector{Float64}}([-1.0], 0.03)

julia> convert(HalfSpace, :(x1 >= -0.03))
HalfSpace{Float64, Vector{Float64}}([-1.0], 0.03)

julia> convert(HalfSpace, :(x1 + x2 <= 2*x4 + 6))
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0, -2.0], 6.0)
```

You can also specify the set of "ambient" variables, even if not
all of them appear:

```jldoctest convert_halfspace
julia> using SymEngine: Basic

julia> convert(HalfSpace, :(x1 + x2 <= 2*x4 + 6), vars=Basic[:x1, :x2, :x3, :x4])
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0, 0.0, -2.0], 6.0)
```
"""
function convert(::Type{HalfSpace{N}}, expr::Expr; vars::Vector{Basic}=Basic[]) where {N}
    @assert _ishalfspace(expr) "the expression $expr does not correspond to a half-space"

    # convert to SymEngine expressions
    linexpr, cmp = _parse_halfspace(expr)

    # check sense of the inequality, assuming < or <= by default (checked before)
    got_geq = cmp ∈ (:(>=), :(>))

    # `a1 x1 + ... + an xn + b [cmp] 0` for [cmp] ∈ {<, <=, >, >=}
    a, b = _parse_linear_expression(linexpr, vars, N)

    return got_geq ? HalfSpace(-a, b) : HalfSpace(a, -b)
end

# type-less default half-space conversion
function convert(::Type{HalfSpace}, expr::Expr; vars::Vector{Basic}=Basic[])
    return convert(HalfSpace{Float64}, expr; vars=vars)
end
