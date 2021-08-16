"""
    require(package::Symbol; fun_name::String="", explanation::String="")

Helper method to check for optional packages and printing an error message.

### Input

- `package`     -- symbol of the package name
- `fun_name`    -- (optional; default: `""`) name of the function that requires
                   the package
- `explanation` -- (optional; default: `""`) additional explanation in the error
                   message

### Output

If the package is loaded, this function has no effect.
Otherwise it prints an error message.

### Algorithm

This function uses `@assert` and hence loses its ability to print an error
message if assertions are deactivated.
"""
function require(package::Symbol; fun_name::String="", explanation::String="")
    @assert isdefined(@__MODULE__, package) "package '$package' not loaded" *
        (fun_name == "" ? "" :
            " (it is required for executing `$fun_name`" *
            (explanation == "" ? "" : " " * explanation) * ")")
end
