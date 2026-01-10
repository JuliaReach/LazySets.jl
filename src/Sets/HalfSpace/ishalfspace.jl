function load_SymEngine_ishalfspace()
    return quote
        using ..LazySets: _is_linear_combination

        """
            _ishalfspace(expr::Expr)

        Determine whether the given expression corresponds to a half-space.

        ### Input

        - `expr` -- a symbolic expression

        ### Output

        `true` if `expr` corresponds to a half-space or `false` otherwise.

        ### Examples

        ```jldoctest
        julia> using LazySets: _ishalfspace

        julia> all(_ishalfspace.([:(x1 <= 0), :(x1 < 0), :(x1 > 0), :(x1 >= 0)]))
        true

        julia> _ishalfspace(:(x1 = 0))
        false

        julia> _ishalfspace(:(2*x1 <= 4))
        true

        julia> _ishalfspace(:(6.1 <= 5.3*f - 0.1*g))
        true

        julia> _ishalfspace(:(2*x1^2 <= 4))
        false

        julia> _ishalfspace(:(x1^2 > 4*x2 - x3))
        false

        julia> _ishalfspace(:(x1 > 4*x2 - x3))
        true
        ```
        """
        function _ishalfspace(expr::Expr)::Bool
            # check that there are three arguments:
            # the comparison symbol, the left-hand side and the right-hand side
            if (length(expr.args) != 3) || !(expr.head == :call)
                return false
            end

            # convert to SymEngine expression
            linexpr, cmp = _parse_halfspace(expr)

            # check that this is an inequality
            if cmp ∉ [:(<=), :(<), :(>=), :(>)]
                return false
            end

            # check if the expression defines a half-space
            return _is_linear_combination(linexpr)
        end
    end
end  # load_SymEngine_ishalfspace

function load_Symbolics_ishalfspace()
    return quote
        if isdefined(Symbolics.SymbolicUtils, :Symbolic)
            using .Symbolics.SymbolicUtils: Symbolic
            const BasicSymbolic2 = Symbolic
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
