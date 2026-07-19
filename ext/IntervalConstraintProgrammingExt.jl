module IntervalConstraintProgrammingExt

using IntervalBoxes: IntervalBox
using IntervalConstraintProgramming: Separator
using Symbolics: @variables  # NOTE: used for metaprogramming below
import LazySets.Approximations: _contract_zonotope_halfspace_ICP

function _contract_zonotope_halfspace_ICP(e::Expr, X::Vector, vars_string::String)
    prefix = "@variables "
    sym = Meta.parse(prefix * vars_string)
    vars = eval(quote
                    $sym
                end)
    separator = eval(quote
                         Separator($e, $vars)
                     end)
    # smallest box containing all points in domain X satisfying constraint
    # (`invokelatest` to avoid world-age issue)
    boundary, _, _ = invokelatest(separator, IntervalBox(X...))  # NOTE: this is an internal function
    return boundary
end

end  # module
