function ishyperplanar(::Hyperplane)
    return true
end

function load_Symbolics_ishyperplanar()
    return quote
        if isdefined(Symbolics.SymbolicUtils, :Symbolic)
            using .Symbolics.SymbolicUtils: Symbolic
            const BasicSymbolic2 = Symbolic
        else
            using .Symbolics: BasicSymbolic
            const BasicSymbolic2 = BasicSymbolic
        end

        # returns `(true, sexpr)` if expr represents a hyperplane,
        # where sexpr is the simplified expression sexpr := LHS - RHS == 0
        # otherwise returns `(false, expr)`
        function _ishyperplanar(expr::BasicSymbolic2)
            got_hyperplane = Symbolics.operation(expr) == ==
            if got_hyperplane
                # simplify to the form a*x + b == 0
                a, b = Symbolics.arguments(expr)
                sexpr = Symbolics.simplify(a - b)
            end
            return got_hyperplane ? (true, sexpr) : (false, expr)
        end
    end
end  # load_Symbolics_ishyperplanar
