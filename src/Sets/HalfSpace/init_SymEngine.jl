import .SymEngine: free_symbols

function free_symbols(expr::Expr, ::Type{<:HalfSpace})
    # get sides of the inequality
    lhs, rhs = convert(SymEngine.Basic, expr.args[2]), convert(SymEngine.Basic, expr.args[3])
    return SymEngine.free_symbols(lhs - rhs)
end

eval(load_SymEngine_ishalfspace())
eval(load_SymEngine_convert_HalfSpace())
