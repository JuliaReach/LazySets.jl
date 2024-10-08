function ishyperplanar(::Hyperplane)
    return true
end

function load_SymEngine_ishyperplanar()
    return quote
        using ..LazySets: _is_linear_combination

        """
            _ishyperplanar(expr::Expr)

        Determine whether the given expression corresponds to a hyperplane.

        ### Input

        - `expr` -- a symbolic expression

        ### Output

        `true` if `expr` corresponds to a half-space or `false` otherwise.

        ### Examples

        ```jldoctest
        julia> using LazySets: _ishyperplanar

        julia> _ishyperplanar(:(x1 = 0))
        true

        julia> _ishyperplanar(:(x1 <= 0))
        false

        julia> _ishyperplanar(:(2*x1 = 4))
        true

        julia> _ishyperplanar(:(6.1 = 5.3*f - 0.1*g))
        true

        julia> _ishyperplanar(:(2*x1^2 = 4))
        false

        julia> _ishyperplanar(:(x1^2 = 4*x2 - x3))
        false

        julia> _ishyperplanar(:(x1 = 4*x2 - x3))
        true
        ```
        """
        function _ishyperplanar(expr::Expr)::Bool
            # check that the head is `=` and there are two arguments:
            # the left-hand side and the right-hand side
            if (length(expr.args) != 2) || !(expr.head == :(=))
                return false
            end

            # convert to SymEngine expression
            linexpr = _parse_hyperplane(expr)

            # check if the expression defines a hyperplane
            return _is_linear_combination(linexpr)
        end
    end
end  # load_SymEngine_ishyperplanar

function load_Symbolics_ishyperplanar()
    return quote
        using .Symbolics: Symbolic

        # returns `(true, sexpr)` if expr represents a hyperplane,
        # where sexpr is the simplified expression sexpr := LHS - RHS == 0
        # otherwise returns `(false, expr)`
        function _ishyperplanar(expr::Symbolic)
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
