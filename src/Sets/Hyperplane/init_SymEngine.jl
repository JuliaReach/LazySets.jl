import .SymEngine: free_symbols

function _parse_hyperplane(expr::Expr)
    lhs = convert(SymEngine.Basic, expr.args[1])
    rhs = :args in fieldnames(typeof(expr.args[2])) ?
          convert(SymEngine.Basic, expr.args[2].args[2]) :
          convert(SymEngine.Basic, expr.args[2])
    return lhs - rhs
end

function free_symbols(expr::Expr, ::Type{<:Hyperplane})
    linexpr = _parse_hyperplane(expr)
    return SymEngine.free_symbols(linexpr)
end

eval(load_SymEngine_ishyperplanar())
eval(load_SymEngine_convert_Hyperplane())
