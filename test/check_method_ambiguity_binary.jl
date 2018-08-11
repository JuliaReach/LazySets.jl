"""
    check_method_ambiguity_binary(op;
                                  [print_results]::Bool=false,
                                  [print_warnings::Bool=print_results])::Bool

Check that a given binary operation does not have method ambiguities for some
combination of set types.

### Input

- `op` -- binary set function name
- `print_results` -- (optional, default: `false`) flag for printing intermediate
                     results
- `print_warnings` -- (optional, default: `print_results`) flag for printing
                      warnings

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
                                       print_warnings::Bool=print_results)::Bool
    types = [AbstractHPolygon{Float64}, AbstractHyperrectangle{Float64},
             AbstractPointSymmetric{Float64},
             AbstractPointSymmetricPolytope{Float64}, AbstractPolygon{Float64},
             AbstractPolytope{Float64}, AbstractSingleton{Float64}]

    polytope_constraints = [LinearConstraint(ones(2), 1.),
                            LinearConstraint([1., -1.], 1.)]
    type2instance = Dict{Type, Vector{LazySet{Float64}}}([
        (AbstractHPolygon{Float64}, [HPolygon{Float64}(polytope_constraints),
                            HPolygonOpt{Float64}(polytope_constraints, 1)]),
        (AbstractHyperrectangle{Float64}, [BallInf(zeros(2), 1.),
                                  Hyperrectangle(zeros(2), ones(2))]),
        (AbstractPointSymmetric{Float64}, [Ball2(zeros(2), 1.),
                                  Ballp(1.5, zeros(2), 1.),
                                  Ellipsoid([1. 0.; 0. 1.])]),
        (AbstractPointSymmetricPolytope{Float64}, [Ball1(zeros(2), 1.),
                                          Zonotope(zeros(2), [1. 0.; 0. 1.])]),
        (AbstractPolygon{Float64}, [VPolygon([zeros(2)])]),
        (AbstractPolytope{Float64}, [HPolytope{Float64}(polytope_constraints)]),
        (AbstractSingleton{Float64}, [Singleton(zeros(2)), ZeroSet(2)])
    ])

    ops = "$op"[10:end] # remove prefix 'LazySets.'
    has_ambiguities = false
    has_other_problems = false
    print_results && println("testing operation $ops")
    for t1 in types
        print_results && println("type 1: $t1")
        for t2 in types
            print_results && println("type 2: $t2")
            if !method_exists(op, (t1, t2))
                if print_warnings
                    warn("no interface operation '$ops(::$(t1), ::$(t2)' exists")
                end
                continue
            end
            print_results && println("interface operation exists")
            for i1 in type2instance[t1]
                print_results && println("instance 1: $(typeof(i1))")
                for i2 in type2instance[t2]
                    print_results && println("instance 2: $(typeof(i2))")
                    try
                        op(i1, i2)
                        print_results && println("test succeeded")
                    catch ex
                        # print error message
#                         showerror(STDOUT, ex); println()

                        # get error message
                        buf = IOBuffer()
                        showerror(buf, ex)
                        msg = String(take!(buf))

                        if startswith(msg, "MethodError: no method matching $ops")
                            # undefined method error, ignore
                            if print_warnings
                                warn("no concrete operation " *
                                     "'$ops(::$(typeof(i1)), " *
                                     "::$(typeof(i2)))' exists")
                            end
                        elseif contains(msg, "is ambiguous")
                            # ambiguity error, remember
                            has_ambiguities = true
                            if print_warnings
                                warn("ambiguous method " *
                                     "'$ops(::$(typeof(i1)), " *
                                     "::$(typeof(i2)))' found; " *
                                     "printing full message")
                                warn(msg)
                            end
                        else
                            # other error, remember
                            has_other_problems = true
                            if print_warnings
                                warn("unknown error message:")
                                warn(msg)
                            end
                        end
                    end
                end
            end
        end
    end
    result = !has_ambiguities
#     result &= !has_other_problems # also report other errors
    return result
end
