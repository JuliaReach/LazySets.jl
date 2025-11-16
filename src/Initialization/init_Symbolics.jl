using .Symbolics: Num  # variable like, e.g. x[1]

# reduce for several variables e.g. when vars = @variables x[1:3] t
_vec(vars::Vector{Any}) = reduce(vcat, vars)
_vec(vars::Vector{Num}) = vars
_vec(vars::Vector{Symbolics.Arr{Num,1}}) = reduce(vcat, vars)
_vec(vars::Vector{Vector{Num}}) = reduce(vcat, vars)
_vec(vars::Vector{Real}) = reduce(vcat, vars)

vSymbolics = pkgversion(Symbolics)
if vSymbolics < v"6.1.0"
    _get_variables(expr::Num) = convert(Vector{Num}, Symbolics.get_variables(expr))
elseif vSymbolics < v"7.0.0"
    # `sort` argument was introduced in Symbolics v6.1
    _get_variables(expr::Num) = convert(Vector{Num}, Symbolics.get_variables(expr; sort=true))
else
    # `sort` argument removed in v7, so manual sorting is needed
    function _get_variables(expr::Num)
        vars_unsorted = Symbolics.get_variables(expr)
        vars_sorted = collect(vars_unsorted)
        sort!(vars_sorted; by=Symbolics.SymbolicUtils.get_degrees)  # NOTE: this is an internal function
        return convert(Vector{Num}, vars_sorted)
    end
end
function _get_variables(expr::Vector{<:Num})
    return unique(reduce(vcat, _get_variables(ex) for ex in expr))
end
