function is_hyperplanar(::Hyperplane)
    return true
end

function load_SymEngine_ishyperplanar()
    return quote
        using .SymEngine: Basic
        using ..LazySets: _is_linearcombination

        """
            _is_hyperplane(expr::Expr)

        Determine whether the given expression corresponds to a hyperplane.

        ### Input

        - `expr` -- a symbolic expression

        ### Output

        `true` if `expr` corresponds to a half-space or `false` otherwise.

        ### Examples

        ```jldoctest
        julia> using LazySets: _is_hyperplane

        julia> _is_hyperplane(:(x1 = 0))
        true

        julia> _is_hyperplane(:(x1 <= 0))
        false

        julia> _is_hyperplane(:(2*x1 = 4))
        true

        julia> _is_hyperplane(:(6.1 = 5.3*f - 0.1*g))
        true

        julia> _is_hyperplane(:(2*x1^2 = 4))
        false

        julia> _is_hyperplane(:(x1^2 = 4*x2 - x3))
        false

        julia> _is_hyperplane(:(x1 = 4*x2 - x3))
        true
        ```
        """
        function _is_hyperplane(expr::Expr)::Bool

            # check that there are three arguments
            # these are the comparison symbol, the left hand side and the right hand side
            if (length(expr.args) != 2) || !(expr.head == :(=))
                return false
            end

            # convert to symengine expressions
            lhs = convert(Basic, expr.args[1])

            if :args in fieldnames(typeof(expr.args[2]))
                # treats the 4 in :(2*x1 = 4)
                rhs = convert(Basic, expr.args[2].args[2])
            else
                rhs = convert(Basic, expr.args[2])
            end

            # check if the expression defines a hyperplane
            return _is_linearcombination(lhs) && _is_linearcombination(rhs)
        end
    end
end  # load_SymEngine_ishyperplanar

function load_Symbolics_ishyperplanar()
    return quote
        using .Symbolics: Symbolic

        # returns `(true, sexpr)` if expr represents a hyperplane,
        # where sexpr is the simplified expression sexpr := LHS - RHS == 0
        # otherwise returns `(false, expr)`
        function _is_hyperplane(expr::Symbolic)
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
