using .SymEngine: Basic, subs

import .SymEngine: free_symbols

"""
   _is_linearcombination(L::Basic)

Determine whether the expression `L` is a linear combination of its symbols.

### Input

- `L` -- expression

### Output

`true` if `L` is a linear combination or false otherwise.

### Examples

```jldoctest
julia> using LazySets: _is_linearcombination

julia> _is_linearcombination(:(2*x1 - 4))
true

julia> _is_linearcombination(:(6.1 - 5.3*f - 0.1*g))
true

julia> _is_linearcombination(:(2*x1^2))
false

julia> _is_linearcombination(:(x1^2 - 4*x2 + x3 + 2))
false
```
"""
function _is_linearcombination(L::Basic)
    return all(isempty.(free_symbols.(diff.(L, free_symbols(L)))))
end

_is_linearcombination(L::Expr) = _is_linearcombination(convert(Basic, L))

"""
    free_symbols(expr::Expr, set_type::Type{LazySet})

Return the free symbols in an expression that represents a given set type.

### Input

- `expr` -- symbolic expression

### Output

A list of symbols, in the form of SymEngine `Basic` objects.

### Examples

```jldoctest
julia> using LazySets: free_symbols

julia> free_symbols(:(x1 <= -0.03), HalfSpace)
1-element Vector{SymEngine.Basic}:
 x1

julia> free_symbols(:(x1 + x2 <= 2*x4 + 6), HalfSpace)
3-element Vector{SymEngine.Basic}:
 x2
 x1
 x4
```
"""
function free_symbols(::Expr, ::Type{<:LazySet}) end

function free_symbols(expr::Expr)
    if _is_hyperplane(expr)
        return free_symbols(expr, HyperPlane)
    elseif _is_halfspace(expr)
        return free_symbols(expr, Halfspace)
    else
        error("the free symbols for the expression $(expr) is not implemented")
    end
end

eval(load_symengine_hyperplane())
eval(load_symengine_halfspace())
