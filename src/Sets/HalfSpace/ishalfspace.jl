function load_SymEngine_ishalfspace()
    return quote
        using .SymEngine: Basic
        import .SymEngine: free_symbols
        using ..LazySets: _is_linearcombination

        """
            _is_halfspace(expr::Expr)

        Determine whether the given expression corresponds to a half-space.

        ### Input

        - `expr` -- a symbolic expression

        ### Output

        `true` if `expr` corresponds to a half-space or `false` otherwise.

        ### Examples

        ```jldoctest
        julia> using LazySets: _is_halfspace

        julia> all(_is_halfspace.([:(x1 <= 0), :(x1 < 0), :(x1 > 0), :(x1 >= 0)]))
        true

        julia> _is_halfspace(:(x1 = 0))
        false

        julia> _is_halfspace(:(2*x1 <= 4))
        true

        julia> _is_halfspace(:(6.1 <= 5.3*f - 0.1*g))
        true

        julia> _is_halfspace(:(2*x1^2 <= 4))
        false

        julia> _is_halfspace(:(x1^2 > 4*x2 - x3))
        false

        julia> _is_halfspace(:(x1 > 4*x2 - x3))
        true
        ```
        """
        function _is_halfspace(expr::Expr)::Bool

            # check that there are three arguments
            # these are the comparison symbol, the left hand side and the right hand side
            if (length(expr.args) != 3) || !(expr.head == :call)
                return false
            end

            # check that this is an inequality
            if !(expr.args[1] in [:(<=), :(<), :(>=), :(>)])
                return false
            end

            # convert to symengine expressions
            lhs, rhs = convert(Basic, expr.args[2]), convert(Basic, expr.args[3])

            # check if the expression defines a half-space
            return _is_linearcombination(lhs) && _is_linearcombination(rhs)
        end
    end
end  # load_SymEngine_ishalfspace

function load_Symbolics_ishalfspace()
    return quote
        using .Symbolics: Symbolic

        # returns `(true, sexpr)` if expr represents a half-space,
        # where sexpr is the simplified expression sexpr := LHS - RHS <= 0
        # otherwise, returns `(false, expr)`
        function _is_halfspace(expr::Symbolic)
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
