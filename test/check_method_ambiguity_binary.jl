"""
    check_method_ambiguity_binary(op;
                                  [print_results]::Bool=false,
                                  [print_warnings::Bool=print_results])::Bool

Check that a given binary operation does not have method ambiguities for some
combination of set types.

### Input

- `op`             -- binary set function name
- `print_results`  -- (optional, default: `false`) flag for printing
                      intermediate results
- `print_warnings` -- (optional, default: `print_results`) flag for printing
                      warnings
- `ignore_errors`  -- (optional, default: `true`) flag for returning `false` if
                      any error occurred

### Output

`true` iff there is no method ambiguity for any set type combination.

### Notes

It is not required that every set type combination is implemented.
The check is only applied to implemented combinations.

This function can also print all intermediate results to STDOUT.
Both method ambiguities/other errors and missing implementations are then
reported as warnings.
Optionally only the warnings can be printed.

### Examples

```julia
julia> check_method_ambiguity_binary(issubset)
```
"""
function check_method_ambiguity_binary(op;
                                       print_results::Bool=false,
                                       print_ambiguity_warnings::Bool=true,
                                       print_warnings::Bool=print_results,
                                       ignore_errors::Bool=true)::Bool
    H1 = LinearConstraint([1.], 1.)
    polytope_constraints_1D = [H1,
                               LinearConstraint([-1.], 0.)]
    sets_1D = LazySet{Float64}[
        Ball1(zeros(1), 1.),
        Ball2(zeros(1), 1.),
        BallInf(zeros(1), 1.),
        Ballp(1.5, zeros(1), 1.),
        Ellipsoid(hcat([1.])),
        EmptySet(),
        H1,
        HPolytope(polytope_constraints_1D),
        HPolyhedron([H1]),
        Hyperplane(ones(1), 1.),
        Hyperrectangle(zeros(1), ones(1)),
        Interval(0., 1.),
        Singleton(zeros(1)),
        VPolytope([zeros(1)]),
        ZeroSet(1),
        Zonotope(zeros(1), hcat([1.; 2.]))
    ]
    H2 = LinearConstraint(ones(2), 1.)
    polytope_constraints_2D = [H2,
                               LinearConstraint([-1., 0.], 0.),
                               LinearConstraint([0., -1.], 0.)]
    sets_2D = LazySet{Float64}[
        Ball1(zeros(2), 1.),
        Ball2(zeros(2), 1.),
        BallInf(zeros(2), 1.),
        Ballp(1.5, zeros(2), 1.),
        Ellipsoid([1. 0.; 0. 1.]),
        EmptySet(),
        H2,
        HPolygon(polytope_constraints_2D),
        HPolygonOpt(polytope_constraints_2D, 1),
        HPolytope(polytope_constraints_2D),
        HPolyhedron([H2]),
        Hyperplane(ones(2), 1.),
        Hyperrectangle(zeros(2), ones(2)),
        Line(ones(2), 1.),
        LineSegment(zeros(2), ones(2)),
        Singleton(zeros(2)),
        Universe(2),
        VPolygon([zeros(2)]),
        VPolytope([zeros(2)]),
        ZeroSet(2),
        Zonotope(zeros(2), [1. 0.; 0. 1.])
    ]

    if startswith("$op", "LazySets.")
        ops = "$op"[10:end] # remove prefix 'LazySets.'
    else
        ops = op
    end
    has_ambiguities = false
    has_other_problems = false
    print_results && println("testing operation $ops")
    for types in [sets_1D, sets_2D]
        for i1 in types
            t1 = typeof(i1)
            print_results && println("type 1: $t1")
            for i2 in types
                t2 = typeof(i2)
                print_results && println("type 2: $t2")
                try
                    op(i1, i2)
                    print_results && println("test succeeded")
                catch ex
                    # print error message
#                     showerror(STDOUT, ex); println()

                    # get error message
                    buf = IOBuffer()
                    showerror(buf, ex)
                    msg = String(take!(buf))

                    if occursin("is ambiguous", msg)
                        # ambiguity error
                        has_ambiguities = true
                        print_warnings && @warn "ambiguous method " *
                            "'$ops(::$(typeof(i1)), ::$(typeof(i2)))' found"
                        print_ambiguity_warnings && @warn msg
                    elseif occursin("no method matching", msg)
                        # not implemented
                        print_results && println("no method " *
                            "'$ops(::$(i1), ::$(i2)' exists")
                    else
                        # other error
                        has_other_problems = true
                        if print_warnings
                            @warn "unknown error message:"
                            @warn msg
                        end
                    end
                end
            end
        end
    end
    result = !has_ambiguities
    if !ignore_errors
        result &= !has_other_problems # also report other errors
    end
    return result
end
