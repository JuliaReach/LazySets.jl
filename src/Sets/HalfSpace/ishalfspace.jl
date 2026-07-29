function load_Symbolics_ishalfspace()
    return quote
        if isdefined(Symbolics.SymbolicUtils, :Symbolic)
            using .Symbolics.SymbolicUtils: Symbolic
            const BasicSymbolic2 = Symbolic
        elseif isdefined(Symbolics.SymbolicUtils, :BasicSymbolic)
            using .Symbolics.SymbolicUtils: BasicSymbolic
            const BasicSymbolic2 = BasicSymbolic
        else
            using .Symbolics: BasicSymbolic
            const BasicSymbolic2 = BasicSymbolic
        end

        # returns `(true, sexpr)` if expr represents a half-space,
        # where sexpr is the simplified expression sexpr := LHS - RHS <= 0
        # otherwise, returns `(false, expr)`
        function _ishalfspace(expr::BasicSymbolic2)
            got_halfspace = true

            # find sense and normalize
            op = Symbolics.operation(expr)
            args = Symbolics.arguments(expr)
            if op in (<=, <)
                a, b = args
                sexpr = Symbolics.simplify(a - b)

            elseif op in (>=, >)
                a, b = args
                sexpr = Symbolics.simplify(b - a)

            elseif (op == |) && (Symbolics.operation(args[1]) == <)
                a, b = Symbolics.arguments(args[1])
                sexpr = Symbolics.simplify(a - b)

            elseif (op == |) && (Symbolics.operation(args[2]) == <)
                a, b = Symbolics.arguments(args[2])
                sexpr = Symbolics.simplify(a - b)

            elseif (op == |) && (Symbolics.operation(args[1]) == >)
                a, b = Symbolics.arguments(args[1])
                sexpr = Symbolics.simplify(b - a)

            elseif (op == |) && (Symbolics.operation(args[2]) == >)
                a, b = Symbolics.arguments(args[2])
                sexpr = Symbolics.simplify(b - a)

            else
                got_halfspace = false
            end

            return got_halfspace ? (true, sexpr) : (false, expr)
        end
    end
end  # load_Symbolics_ishalfspace
