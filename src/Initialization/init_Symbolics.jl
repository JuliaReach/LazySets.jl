using .Symbolics: Num  # variable like, e.g. x[1]

# reduce for several variables e.g. when vars = @variables x[1:3] t
_vec(vars::Vector{Any}) = reduce(vcat, vars)
_vec(vars::Vector{Num}) = vars
_vec(vars::Vector{Symbolics.Arr{Num,1}}) = reduce(vcat, vars)
_vec(vars::Vector{Vector{Num}}) = reduce(vcat, vars)
_vec(vars::Vector{Real}) = reduce(vcat, vars)

# `sort` argument was introduced in Symbolics v6.1
@static if VERSION < v"1.10"
    _get_variables(expr::Num) = convert(Vector{Num}, Symbolics.get_variables(expr))
else
    _get_variables(expr::Num) = convert(Vector{Num}, Symbolics.get_variables(expr; sort=true))
end
function _get_variables(expr::Vector{<:Num})
    return unique(reduce(vcat, _get_variables(ex) for ex in expr))
end
