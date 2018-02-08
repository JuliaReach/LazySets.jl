export sign_cadlag,
       check_method_implementation,
       check_method_ambiguity_binary,
       @commutative_neutral,
       @commutative_absorbing

"""
    sign_cadlag(x::N)::N where {N<:Real}

This function works like the sign function but is ``1`` for input ``0``.

### Input

- `x` -- real scalar

### Output

``1`` if ``x ≥ 0``, ``-1`` otherwise.

### Notes

This is the sign function right-continuous at zero (see
[càdlàg function](https://en.wikipedia.org/wiki/C%C3%A0dl%C3%A0g)).
It can be used with vector-valued arguments via the dot operator.

### Examples

```jldoctest
julia> sign_cadlag.([-0.6, 1.3, 0.0])
3-element Array{Float64,1}:
 -1.0
  1.0
  1.0
```
"""
function sign_cadlag(x::N)::N where {N<:Real}
    return x < zero(x) ? -one(x) : one(x)
end

"""
    check_method_implementation(interface::Type,
                                func_name,
                                args_funcs::AbstractVector{Function},
                                [print_results]::Bool=false
                               )::Bool

Check that a given (interface) function is implemented by all subtypes.

### Input

- `interface`     -- parent interface type
- `func_name`     -- function name
- `args_funcs`    -- list of functions that each map a type to an argument
                     signature tuple
- `print_results` -- (optional, default: `false`) flag for printing intermediate
                     results

### Output

`true` iff all subtypes implement the given function.

### Notes

It is sufficient that a subinterface implements the function, i.e., it is not
required that every *concrete* type implements the function.

This function can also print all intermediate results to STDOUT.

### Examples

```jldoctest
julia> check_method_implementation(LazySet, σ,
        Function[S -> (AbstractVector{Float64}, S)])
true
```
"""
function check_method_implementation(interface::Type,
                                     func_name,
                                     args_funcs::AbstractVector{Function};
                                     print_results::Bool=false
                                    )::Bool
    has_subtypes = false # NOTE: 'isleaftype' does not work (type parameters)
    for subtype in subtypes(interface)
        has_subtypes = true
        found = false
        for args_func in args_funcs
            if method_exists(func_name, args_func(subtype))
                if print_results
                    println("found implementation of $func_name for $subtype")
                end
                found = true
                break
            end
        end
        if !found && !check_method_implementation(subtype, func_name,
                                                  args_funcs,
                                                  print_results=print_results)
            if print_results
                println("no implementation of $func_name for $subtype")
            end
            return false
        end
    end
    return has_subtypes
end

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
julia> check_method_ambiguity_binary(is_subset)
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
                                  Ellipsoid(eye(2))]),
        (AbstractPointSymmetricPolytope{Float64}, [Ball1(zeros(2), 1.),
                                          Zonotope(zeros(2), eye(2))]),
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

"""
    @commutative_neutral(SET, NEUT)

Creates functions to make a set type behave commutative with a given neutral
element set type.

### Input

- `SET`  -- set type
- `NEUT` -- set type of the neutral element

### Output

Three function definitions.

### Examples

`@commutative_neutral(ConvexHull, EmptySet)` creates the following functions:
* `ConvexHull(X::LazySet{N}, ::EmptySet{N}) where {N<:Real} = X`
* `ConvexHull(::EmptySet{N}, X::LazySet{N}) where {N<:Real} = X`
* `ConvexHull(Y::EmptySet{N}, ::EmptySet{N}) where {N<:Real} = Y`
"""
macro commutative_neutral(SET, NEUT)
    @eval begin
        function $SET(X::LazySet{N}, ::$NEUT{N}) where {N<:Real}
            return X
        end
        function $SET(::$NEUT{N}, X::LazySet{N}) where {N<:Real}
            return X
        end
        function $SET(Y::$NEUT{N}, ::$NEUT{N}) where {N<:Real}
            return Y
        end
    end
    return nothing
end

"""
    @commutative_absorbing(SET, ABS)

Creates functions to make a set type behave commutative with a given absorbing
element set type.

### Input

- `SET`  -- set type
- `ABS` -- set type of the absorbing element

### Output

Three function definitions.

### Examples

`@commutative_absorbing(MinkowskiSum, EmptySet)` creates the following
functions:
* `MinkowskiSum(::LazySet{N}, Y::EmptySet{N}) where {N<:Real} = Y`
* `MinkowskiSum(Y::EmptySet{N}, ::LazySet{N}) where {N<:Real} = Y`
* `MinkowskiSum(Y::EmptySet{N}, ::EmptySet{N}) where {N<:Real} = Y`
"""
macro commutative_absorbing(SET, ABS)
    @eval begin
        function $SET(::LazySet{N}, Y::$ABS{N}) where {N<:Real}
            return Y
        end
        function $SET(Y::$ABS{N}, ::LazySet{N}) where {N<:Real}
            return Y
        end
        function $SET(Y::$ABS{N}, ::$ABS{N}) where {N<:Real}
            return Y
        end
    end
    return nothing
end
