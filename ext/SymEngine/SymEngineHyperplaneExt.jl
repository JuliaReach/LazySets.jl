using LazySets.HyperplaneModule: Hyperplane
using SymEngine: Basic
import Base: convert
import SymEngine: free_symbols

function _parse_hyperplane(expr::Expr)
    lhs = convert(Basic, expr.args[1])
    rhs = :args in fieldnames(typeof(expr.args[2])) ?
          convert(Basic, expr.args[2].args[2]) :
          convert(Basic, expr.args[2])
    return lhs - rhs
end

function free_symbols(expr::Expr, ::Type{<:Hyperplane})
    linexpr = _parse_hyperplane(expr)
    return SymEngine.free_symbols(linexpr)
end

"""
    _ishyperplanar(expr::Expr)

Determine whether the given expression corresponds to a hyperplane.

### Input

- `expr` -- a symbolic expression

### Output

`true` if `expr` corresponds to a half-space or `false` otherwise.

### Examples

```jldoctest
julia> using LazySets

julia> import SymEngine

julia> SymEngineExt = Base.get_extension(LazySets, :SymEngineExt);

julia> using .SymEngineExt: _ishyperplanar

julia> _ishyperplanar(:(x1 = 0))
true

julia> _ishyperplanar(:(x1 <= 0))
false

julia> _ishyperplanar(:(2*x1 = 4))
true

julia> _ishyperplanar(:(6.1 = 5.3*f - 0.1*g))
true

julia> _ishyperplanar(:(2*x1^2 = 4))
false

julia> _ishyperplanar(:(x1^2 = 4*x2 - x3))
false

julia> _ishyperplanar(:(x1 = 4*x2 - x3))
true
```
"""
function _ishyperplanar(expr::Expr)::Bool
    # check that the head is `=` and there are two arguments:
    # the left-hand side and the right-hand side
    if (length(expr.args) != 2) || !(expr.head == :(=))
        return false
    end

    # convert to SymEngine expression
    linexpr = _parse_hyperplane(expr)

    # check if the expression defines a hyperplane
    return _is_linear_combination(linexpr)
end

"""
    convert(::Type{Hyperplane{N}}, expr::Expr; vars=Vector{Basic}=Basic[]) where {N}

Return a `LazySet.Hyperplane` given a symbolic expression that represents a hyperplane.

### Input

- `expr` -- a symbolic expression
- `vars` -- (optional, default: `Basic[]`): set of variables with respect to which the
            gradient is taken; if empty, we take the free symbols in the given expression

### Output

A `Hyperplane`, in the form `ax = b`.

### Examples

```jldoctest convert_hyperplane
julia> using LazySets

julia> import SymEngine

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
    @assert _ishyperplanar(expr) "the expression $expr does not correspond to a Hyperplane"

    # convert to SymEngine expression
    linexpr = _parse_hyperplane(expr)

    # a1 x1 + ... + an xn + b = 0
    a, b = _parse_linear_expression(linexpr, vars, N)

    return Hyperplane(a, -b)
end

# type-less default Hyperplane conversion
function convert(::Type{Hyperplane}, expr::Expr; vars::Vector{Basic}=Basic[])
    return convert(Hyperplane{Float64}, expr; vars=vars)
end
