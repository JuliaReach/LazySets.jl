import .SymEngine: free_symbols

function free_symbols(expr::Expr, ::Type{<:Hyperplane})
    # get sides of the equality
    lhs = convert(SymEngine.Basic, expr.args[1])

    # treats the 4 in :(2*x1 = 4)
    rhs = :args in fieldnames(typeof(expr.args[2])) ? convert(SymEngine.Basic, expr.args[2].args[2]) :
            convert(SymEngine.Basic, expr.args[2])
    return free_symbols(lhs - rhs)
end

eval(load_SymEngine_ishyperplanar())
eval(load_SymEngine_convert_Hyperplane())
