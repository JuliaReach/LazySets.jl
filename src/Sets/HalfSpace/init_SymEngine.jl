import .SymEngine: free_symbols

function _parse_halfspace(expr::Expr)
    lhs = convert(SymEngine.Basic, expr.args[2])
    rhs = convert(SymEngine.Basic, expr.args[3])
    cmp = expr.args[1]
    return (lhs - rhs, cmp)
end

function free_symbols(expr::Expr, ::Type{<:HalfSpace})
    linexpr, _ = _parse_halfspace(expr)
    return SymEngine.free_symbols(linexpr)
end

eval(load_SymEngine_ishalfspace())
eval(load_SymEngine_convert_HalfSpace())
