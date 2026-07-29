module LazySetsSymbolicsExt

import Symbolics
using Symbolics: Arr, Num, get_variables
import SymbolicUtils

include("Symbolics/HalfSpaceExt.jl")
include("Symbolics/HPolyhedronExt.jl")
include("Symbolics/HPolytopeExt.jl")
include("Symbolics/HyperplaneExt.jl")

"""
    _vec(vars)

Transform a tuple of operations into one vector of operations.

### Input

- `vars` -- tuple where each element is either variable-like (`Num`) or a
            vector of variables (`Vector{Num}`)

### Output

A vector of `Operation` obtained by concatenating each tuple component.

## Examples

```jldoctest
julia> using LazySets, Symbolics

julia> vars = @variables x[1:2] y
2-element Vector{Any}:
  x[1:2]
 y

julia> LazySetsSymbolicsExt = Base.get_extension(LazySets, :LazySetsSymbolicsExt);

julia> LazySetsSymbolicsExt._vec(vars)
3-element Vector{Num}:
 x[1]
 x[2]
    y
```
"""
function _vec end

# reduce for several variables e.g. when vars = @variables x[1:3] t
_vec(vars::Vector{Any}) = reduce(vcat, vars)
_vec(vars::Vector{Num}) = vars
_vec(vars::Vector{Arr{Num,1}}) = reduce(vcat, vars)
_vec(vars::Vector{Vector{Num}}) = reduce(vcat, vars)
_vec(vars::Vector{Real}) = reduce(vcat, vars)

vSymbolics = pkgversion(Symbolics)
if vSymbolics < v"6.1.0"
    function _get_variables(expr::Num)
        return convert(Vector{Num}, get_variables(expr))
    end
elseif vSymbolics < v"7.0.0"
    # `sort` argument was introduced in Symbolics v6.1
    function _get_variables(expr::Num)
        return convert(Vector{Num}, get_variables(expr; sort=true))
    end
else
    # `sort` argument removed in v7, so manual sorting is needed
    function _get_variables(expr::Num)
        vars_unsorted = get_variables(expr)
        vars_sorted = collect(vars_unsorted)
        sort!(vars_sorted; by=SymbolicUtils.get_degrees)  # NOTE: this is an internal function
        return convert(Vector{Num}, vars_sorted)
    end
end

function _get_variables(expr::Vector{<:Num})
    return unique(reduce(vcat, _get_variables(ex) for ex in expr))
end

end  # module
