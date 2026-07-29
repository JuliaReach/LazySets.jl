using Symbolics: Num
import LazySets.HPolytopeModule: HPolytope

"""
    HPolytope(expr::Vector{<:Num}, vars=_get_variables(expr);
                [N]::Type{<:Real}=Float64, [check_boundedness]::Bool=false)

Return the polytope in constraint representation given by a list of symbolic
expressions.

### Input

- `expr` -- vector of symbolic expressions that describes each constraint
- `vars` -- (optional, default: `_get_variables(expr)`) if an array of variables
            is given, use those as the ambient variables in the set with respect
            to which derivations take place; otherwise, use only the variables
            that appear in the given expression (but be careful because the
            order may be incorrect; it is advised to always pass `vars`
            explicitly)
- `N`    -- (optional, default: `Float64`) the numeric type of the returned
            polytope
- `check_boundedness` -- (optional, default: `false`) flag to check boundedness

### Output

An `HPolytope`.

### Examples

```jldoctest
julia> using LazySets, Symbolics

julia> vars = @variables x y
2-element Vector{Num}:
 x
 y

julia> HPolytope([x <= 1, x >= 0, y <= 1, y >= 0], vars)
HPolytope{Float64, Vector{Float64}}(HalfSpace{Float64, Vector{Float64}}[HalfSpace{Float64, Vector{Float64}}([1.0, 0.0], 1.0), HalfSpace{Float64, Vector{Float64}}([-1.0, 0.0], 0.0), HalfSpace{Float64, Vector{Float64}}([0.0, 1.0], 1.0), HalfSpace{Float64, Vector{Float64}}([0.0, -1.0], 0.0)])
```
"""
function HPolytope(expr::Vector{<:Num}, vars::AbstractVector{Num}=_get_variables(expr);
                   N::Type{<:Real}=Float64, check_boundedness::Bool=false)
    return HPolytope([HalfSpace(ex, vars; N=N) for ex in expr];
                     check_boundedness=check_boundedness)
end

function HPolytope(expr::Vector{<:Num}, vars; N::Type{<:Real}=Float64,
                   check_boundedness::Bool=false)
    return HPolytope(expr, _vec(vars); N=N, check_boundedness=check_boundedness)
end
