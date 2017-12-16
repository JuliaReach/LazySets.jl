import Base.<=

export sign_cadlag,
       jump2pi,
       check_method_implementation

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
    jump2pi(x::Float64)::Float64

Return ``x + 2π`` if ``x`` is negative, otherwise return ``x``.

### Input

- `x` -- real scalar

### Output

``x + 2π`` if ``x`` is negative, ``x`` otherwise.

### Examples

```jldoctest
julia> jump2pi(0.0)
0.0
julia> jump2pi(-0.5)
5.783185307179586
julia> jump2pi(0.5)
0.5
```
"""
function jump2pi(x::Float64)::Float64
    x < 0.0 ? 2.0 * pi + x : x
end

"""
    <=(u::AbstractVector{Float64}, v::AbstractVector{Float64})::Bool

Compares two 2D vectors by their direction.

### Input

- `u` --  first 2D direction
- `v` --  second 2D direction

### Output

True iff ``\\arg(u) [2π] ≤ \\arg(v) [2π]``

### Notes

The argument is measured in counter-clockwise fashion, with the 0 being the
direction (1, 0).
"""
function <=(u::AbstractVector{Float64}, v::AbstractVector{Float64})::Bool
    return jump2pi(atan2(u[2], u[1])) <= jump2pi(atan2(v[2], v[1]))
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

```julia
check_method_implementation(LazySet, σ,
        Function[S -> (AbstractVector{Float64}, S)])
true
```
"""
function check_method_implementation(interface::Type,
                                     func_name,
                                     args_funcs::AbstractVector{Function};
                                     print_results::Bool=false
                                    )::Bool
    for subtype in subtypes(interface)
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
    return true
end
