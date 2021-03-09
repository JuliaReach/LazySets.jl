"""
    check_method_implementation(interface::Type,
                                func_name,
                                args_funcs::AbstractVector{Function};
                                [ignore_types]::Vector{Type}=Type[],
                                [print_results]::Bool=false
                               )::Bool

Check that a given (interface) function is implemented by all subtypes.

### Input

- `interface`     -- parent interface type
- `func_name`     -- function name
- `args_funcs`    -- list of functions that each map a type to an argument
                     signature tuple
- `ignore_types`  -- (optional, default: `Type[]`) list of types that should be
                     ignored
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
                                     ignore_types::Vector{Type}=Type[],
                                     print_results::Bool=false
                                    )::Bool
    # first collect all base types that are subtypes of this interface
    # NOTE: 'isleaftype' does not work due to type parameters
    base_types = LazySets.subtypes(interface, true)

    # now check all base types
    for subtype in base_types
        if subtype ∈ ignore_types
            if print_results
                println("ignoring type $subtype")
            end
            continue
        end
        found = false
        for args_func in args_funcs
            if hasmethod(func_name, args_func(subtype))
                if print_results
                    println("found implementation of $func_name for $subtype")
                end
                found = true
                break
            end
        end
        if !found
            if print_results
                println("no implementation of $func_name for $subtype")
            end
            return false
        end
    end
    return true
end
