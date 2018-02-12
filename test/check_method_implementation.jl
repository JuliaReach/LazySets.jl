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
