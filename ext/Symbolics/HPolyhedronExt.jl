using Symbolics: Num, gradient, substitute, value
import LazySets.HPolyhedronModule: HPolyhedron

"""
    HPolyhedron(expr::Vector{<:Num}, vars=_get_variables(expr);
                N::Type{<:Real}=Float64)

Return a polyhedron in constraint representation given by a list of symbolic
expressions.

### Input

- `expr` -- vector of symbolic expressions that describes each half-space
- `vars` -- (optional, default: `_get_variables(expr)`), if an array of
            variables is given, use those as the ambient variables in the set
            with respect to which derivations take place; otherwise, use only
            the variables that appear in the given expression (but be careful
            because the order may be incorrect; it is advised to always pass
            `vars` explicitly)
- `N`    -- (optional, default: `Float64`) the numeric type of the returned set

### Output

An `HPolyhedron`.

### Examples

```jldoctest
julia> using LazySets, Symbolics

julia> vars = @variables x y
2-element Vector{Num}:
 x
 y

julia> HPolyhedron([x + y <= 1, x + y >= -1], vars)
HPolyhedron{Float64, Vector{Float64}}(HalfSpace{Float64, Vector{Float64}}[HalfSpace{Float64, Vector{Float64}}([1.0, 1.0], 1.0), HalfSpace{Float64, Vector{Float64}}([-1.0, -1.0], 1.0)])

julia> X = HPolyhedron([x == 0, y <= 0], vars)
HPolyhedron{Float64, Vector{Float64}}(HalfSpace{Float64, Vector{Float64}}[HalfSpace{Float64, Vector{Float64}}([1.0, 0.0], -0.0), HalfSpace{Float64, Vector{Float64}}([-1.0, -0.0], 0.0), HalfSpace{Float64, Vector{Float64}}([0.0, 1.0], -0.0)])
```
"""
function HPolyhedron(expr::Vector{<:Num}, vars::AbstractVector{Num}; N::Type{<:Real}=Float64)
    # TODO this method allows hyperplanes but the HPolytope method does not
    clist = Vector{HalfSpace{N,Vector{N}}}()
    sizehint!(clist, length(expr))
    zeroed_vars = Dict(v => zero(N) for v in vars)
    vars_list = collect(vars)
    for ex in expr
        exval = value(ex)
        got_hyperplane, sexpr = _ishyperplanar(exval)
        if !got_hyperplane
            got_halfspace, sexpr = _ishalfspace(exval)
            if !got_halfspace
                throw(ArgumentError("expected an expression describing either " *
                                    "a half-space of a hyperplane, got $expr"))
            end
        end

        coeffs = [N(value(α)) for α in gradient(sexpr, vars_list)]
        β = -N(value(substitute(sexpr, zeroed_vars)))

        push!(clist, HalfSpace(coeffs, β))
        if got_hyperplane
            push!(clist, HalfSpace(-coeffs, -β))
        end
    end
    return HPolyhedron(clist)
end

function HPolyhedron(expr::Vector{<:Num}; N::Type{<:Real}=Float64)
    return HPolyhedron(expr, _get_variables(expr); N=N)
end
function HPolyhedron(expr::Vector{<:Num}, vars; N::Type{<:Real}=Float64)
    return HPolyhedron(expr, _vec(vars); N=N)
end
