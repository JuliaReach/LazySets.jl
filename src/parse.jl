using SymEngine

import Base.convert
import SymEngine.free_symbols

"""
   is_linearcombination(L::Expr)

Return whether the `L` is a linear combination of its symbols.

### Input

- `L` -- symbolic expression

### Output

`true` if `L` is a linear combination and `false` otherwise.

### Examples

```jldoctest
julia> import LazySets.is_linearcombination

julia> is_linearcombination(:(2*x1 - 4))
true

julia> is_linearcombination(:(6.1 - 5.3*f - 0.1*g))
true

julia> is_linearcombination(:(2*x1^2))
false

julia> is_linearcombination(:(x1^2 - 4*x2 + x3 + 2))
false
```
"""
is_linearcombination(L::Expr) = is_linearcombination(convert(Basic, L))

function is_linearcombination(L::Basic)
    return all(isempty.(free_symbols.(diff.(L, free_symbols(L)))))
end

"""
    is_halfspace(expr::Expr)

Return whether the given expression corresponds to a halfspace.

### Input

- `expr` -- symbolic expression

### Output

`true` if `expr` corresponds to a halfspace and `false` otherwise.

### Examples

Note that the half-space doesn't need to be expressed in "canonical form"
`ax ≤ b`:

```jldoctest
julia> import LazySets.is_halfspace

julia> all(is_halfspace.([:(x1 <= 0), :(x1 < 0), :(x1 > 0), :(x1 >= 0)]))
true

julia> is_halfspace(:(2*x1 <= 4))
true

julia> is_halfspace(:(2*x1 ≤ 4))
true

julia> is_halfspace(:(6.1 <= 5.3*f - 0.1*g))
true

julia> is_halfspace(:(2*x1^2 <= 4))
false

julia> is_halfspace(:(x1^2 > 4*x2 - x3))
false

julia> is_halfspace(:(x1 > 4*x2 - x3))
true
```
"""
function is_halfspace(expr::Expr)

    # check that there are three arguments
    # they are: comparison symbol, left hand side and right hand side
    if (length(expr.args) != 3) || !(expr.head == :call)
        return false
    end

    # check that this is an inequality
    if !(expr.args[1] in [:(<=), :(≤), :(<), :(>=), :(≥), :(>)])
        return false
    end

    # convert to symengine expressions
    lhs, rhs = convert(Basic, expr.args[2]), convert(Basic, expr.args[3])

    # check if the expression defines a halfspace
    return is_linearcombination(lhs) && is_linearcombination(rhs)
end

"""
    is_hyperplane(expr::Expr)

Return whether the given expression corresponds to a hyperplane.

### Input

- `expr` -- symbolic expression

### Output

`true` if `expr` corresponds to a halfspace or `false` otherwise.

### Examples

```jldoctest
julia> import LazySets.is_hyperplane

julia> is_hyperplane(:(x1 = 0))
true

julia> is_hyperplane(:(2*x1 = 4))
true

julia> is_hyperplane(:(6.1 = 5.3*f - 0.1*g))
true

julia> is_hyperplane(:(2*x1^2 = 4))
false

julia> is_hyperplane(:(x1^2 = 4*x2 - x3))
false

julia> is_hyperplane(:(x1 = 4*x2 - x3))
true
```
"""
function is_hyperplane(expr::Expr)::Bool

    # check that there are three arguments
    # they are: comparison symbol, left hand side and right hand side
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
    return is_linearcombination(lhs) && is_linearcombination(rhs)
end

"""
    convert(::Type{HalfSpace{N}}, expr::Expr;
            vars::Vector{Symbol}=Symbol[]) where {N}

Return a `HalfSpace` given a symbolic expression that represents a halfspace.

### Input

- `expr` -- symbolic expression
- `vars` -- (optional, default: `nothing`): set of variables with respect to which
            the gradient is taken; if nothing, the free symbols in the given
            expression are considered

### Output

A `HalfSpace` in the form `ax <= b`.

### Examples

```jldoctest convert_halfspace
julia> import LazySets.convert

julia> convert(HalfSpace, :(x1 <= -0.03))
HalfSpace{Float64}([1.0], -0.03)

julia> convert(HalfSpace, :(x1 < -0.03))
HalfSpace{Float64}([1.0], -0.03)

julia> convert(HalfSpace, :(x1 > -0.03))
HalfSpace{Float64}([-1.0], 0.03)

julia> convert(HalfSpace, :(x1 >= -0.03))
HalfSpace{Float64}([-1.0], 0.03)

julia> convert(HalfSpace, :(x1 + x2 <= 2*x4 + 6))
HalfSpace{Float64}([1.0, 1.0, -2.0], 6.0)
```

The set of "ambient" variables can be specified by passing the `vars` argument,
even if not all of them appear in the halfspace's equation:

```jldoctest convert_halfspace
julia> convert(HalfSpace, :(x1 + x2 <= 2*x4 + 6), vars=[:x1, :x2, :x3, :x4])
HalfSpace{Float64}([1.0, 1.0, 0.0, -2.0], 6.0)
```
"""
function convert(::Type{HalfSpace{N}}, expr::Expr; vars::Vector{Symbol}=Symbol[]) where N

    @assert is_halfspace(expr) "the expression $expr does not correspond to a halfspace"

    vars = convert(Vector{Basic}, vars)

    # check sense of the inequality, assuming < or <= by default
    got_geq = expr.args[1] in [:(>=), :(>), :(≥)]

    # get sides of the inequality
    lhs, rhs = convert(Basic, expr.args[2]), convert(Basic, expr.args[3])

    # a1 x1 + ... + an xn + K [cmp] 0 for cmp in <, <=, >, >=
    eq = lhs - rhs
    if isempty(vars)
        vars = free_symbols(eq)
    end
    K = subs(eq, [vi => zero(N) for vi in vars]...)
    a = convert(Basic, eq - K)

    # convert to numeric types
    K = convert(N, K)
    a = convert(Vector{N}, diff.(a, vars))

    if got_geq
        return HalfSpace(-a, K)
    else
        return HalfSpace(a, -K)
    end
end

# type-less default halfspace conversion
# TODO: make it work using dispatch on the coeffs of the expr, see SymEngine#155
convert(::Type{HalfSpace}, expr::Expr; vars::Vector{Symbol}=Symbol[]) = convert(HalfSpace{Float64}, expr; vars=vars)

"""
    convert(::Type{Hyperplane{N}}, expr::Expr;
            vars::Vector{Symbol}=Symbol[]) where {N}

Return a `Hyperplane` given a symbolic expression that represents a hyperplane.

### Input

- `expr` -- symbolic expression
- `vars` -- (optional, default: `nothing`): set of variables with respect to which
            the gradient is taken; if nothing, the free symbols in the given
            expression are considered

### Output

A `Hyperplane`, in the form `ax = b`.

### Examples

```jldoctest convert_hyperplane
julia> import LazySets.convert

julia> convert(Hyperplane, :(x1 = -0.03))
Hyperplane{Float64}([1.0], -0.03)

julia> convert(Hyperplane, :(x1 + 0.03 = 0))
Hyperplane{Float64}([1.0], -0.03)

julia> convert(Hyperplane, :(x1 + x2 = 2*x4 + 6))
Hyperplane{Float64}([1.0, 1.0, -2.0], 6.0)
```

Use the `vars` option to specify the set of "ambient" variables in the
hyperplane:

```jldoctest convert_hyperplane
julia> convert(Hyperplane, :(x1 + x2 = 2*x4 + 6), vars=[:x1, :x2, :x3, :x4])
Hyperplane{Float64}([1.0, 1.0, 0.0, -2.0], 6.0)
```
"""
function convert(::Type{Hyperplane{N}}, expr::Expr; vars::Vector{Symbol}=Symbol[]) where N

    vars = convert(Vector{Basic}, vars)

    @assert is_hyperplane(expr) "the expression $expr does not correspond to a Hyperplane"

    # get sides of the inequality
    lhs = convert(Basic, expr.args[1])

    # treats the 4 in :(2*x1 = 4)
    rhs = :args in fieldnames(typeof(expr.args[2])) ? convert(Basic, expr.args[2].args[2]) : convert(Basic, expr.args[2])

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
convert(::Type{Hyperplane}, expr::Expr; vars::Vector{Symbol}=Symbol[]) = convert(Hyperplane{Float64}, expr; vars=vars)

"""
    free_symbols(expr::Expr, set_type::Type{LazySet})

Return the free symbols in an expression that represents a given set type.

### Input

- `expr` -- symbolic expression

### Output

A list of symbols, in the form of SymEngine `Basic` objects.

### Examples

```jldoctest free_symbols
julia> import SX.free_symbols

julia> free_symbols(:(x1 <= -0.03), HalfSpace)
1-element Array{SymEngine.Basic,1}:
 x1

julia> free_symbols(:(x1 + x2 <= 2*x4 + 6), HalfSpace)
3-element Array{SymEngine.Basic,1}:
 x2
 x1
 x4
```
"""
function free_symbols(expr::Expr, set_type::Type{<:HalfSpace})
    # get sides of the inequality
    lhs, rhs = convert(Basic, expr.args[2]), convert(Basic, expr.args[3])
    return free_symbols(lhs - rhs)
end

function free_symbols(expr::Expr, set_type::Type{<:Hyperplane})
    # get sides of the inequality
    lhs = convert(Basic, expr.args[1])

    # treats the 4 in :(2*x1 = 4)
    rhs = :args in fieldnames(typeof(expr.args[2])) ? convert(Basic, expr.args[2].args[2]) : convert(Basic, expr.args[2])
    return free_symbols(lhs - rhs)
end

function free_symbols(expr::Expr)
    if is_hyperplane(expr)
        return free_symbols(expr, HyperPlane)
    elseif is_halfspace(expr)
        return free_symbols(expr, Halfspace)
    else
        error("the free symbols for the expression $(expr) is not implemented")
    end
end

"""
    check_malformed_expression(s; assignment=false)

Returns the list of expressions corresponding to a given SX string.

### Input

- `s`          -- string
- `assignment` -- (optional, default: `false`)

### Output

Vector of expressions, equations or inequalities.

### Examples

```jldoctest parse_sxmath
julia> import SX.parse_sxmath

julia> parse_sxmath("x >= 0")
1-element Array{Expr,1}:
 :(x >= 0)

julia> parse_sxmath("x' == x & v' == -0.75*v")
2-element Array{Expr,1}:
 :(x' = x)
 :(v' = -0.75v)

julia> parse_sxmath("x == 0 & v <= 0")
2-element Array{Expr,1}:
 :(x = 0)
 :(v <= 0)
```

Parentheses are ignored:

```jldoctest parse_sxmath
julia> parse_sxmath("(x == 0) & (v <= 0)")
2-element Array{Expr,1}:
 :(x = 0)
 :(v <= 0)
```

Splitting is also performend over double ampersand symbols:

```jldoctest parse_sxmath
julia> parse_sxmath("x == 0 && v <= 0")
2-element Array{Expr,1}:
 :(x = 0)
 :(v <= 0)
```

If you want to parse an assignment, use the `assignment` flag:

```jldoctest parse_sxmath
julia> parse_sxmath("x := -x*0.1", assignment=true)
1-element Array{Expr,1}:
 :(x = -x * 0.1)
```

 Check that we can parse expressions involving parentheses:

```jldoctest parse_sxmath
julia> parse_sxmath("(t <= 125 & y>= -100)")
2-element Array{Expr,1}:
 :(t <= 125)
 :(y >= -100)
julia> parse_sxmath("t <= 125 & (y>= -100)")
2-element Array{Expr,1}:
 :(t <= 125)
 :(y >= -100)
```

### Algorithm

First a sanity check (assertion) is made that the expression makes a coherent use
of parentheses.

Then, the following steps are done (in the given order):

- split the string with the `&` key, or `&&`
- remove trailing whitespace of each substring
- replace double `==` with single `=`
- detect unbalanced parentheses (beginning and final subexpressions) and in that case delete them
- cast to a Julia expression with `parse`

### Notes

For assignments, the nomenclature `:=` is also valid and here it is replaced
to `=`, but you need to set `assignment=true` for this replacement to take effect.

The characters `'('` and `')'` are deleted (replaced by the empty character),
whenever it is found that there are unbalanced parentheses after the expression is
splitted into subexpressions.
"""
function check_malformed_expression(s; assignment=false)
    accept = true

    count_left_parentheses = s -> count(c -> c == '(', collect(s))
    count_right_parentheses = s -> count(c -> c == ')', collect(s))
    if count_left_parentheses(s) != count_right_parentheses(s)
        return (false, "the expression $(s) is not well formed; parentheses do not match")
    end

    # canonicalize expression
    expr = strip.(split(replace(s, "&&" => "&"), "&"))
    if assignment
        expr = replace.(expr, Ref(":=" => "="))
    end
    expr = replace.(expr, Ref("==" => "="))

    # remove irrelevant parentheses from the beginning and the end
    for (i, expri) in enumerate(expr)
        m = count_left_parentheses(expri)
        n = count_right_parentheses(expri)
        if m == n
            nothing
        elseif m > n && expri[1] == '('
            expr[i] = expri[2:end]
        elseif m < n && expri[end] == ')'
            expr[i] = expri[1:end-1]
        else
            error("malformed subexpression $(expri)")
        end
    end

    return Base.Meta.parse.(expr)
end