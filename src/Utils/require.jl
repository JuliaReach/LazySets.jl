"""
    require(package::Symbol; [fun_name]::String="", [explanation]::String="")

Helper method to check for optional packages and print an error message.

### Input

- `package`     -- symbol (the package name)
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
    check = isdefined(@__MODULE__, package)
    @assert check "package '$package' not loaded" *
        (fun_name == "" ? "" :
            " (it is required for executing `$fun_name`" *
            (explanation == "" ? "" : " " * explanation) * ")")
end

"""
    require(packages::Vector{Symbol}; [fun_name]::String="",
                                      [explanation]::String="")

Helper method to check for optional packages and print an error message.

### Input

- `packages`    -- list of symbols (the package names)
- `fun_name`    -- (optional; default: `""`) name of the function that requires
                   the package
- `explanation` -- (optional; default: `""`) additional explanation in the error
                   message

### Output

If at least one of the packages is loaded, this function has no effect.
Otherwise it prints an error message.

### Algorithm

This function uses `@assert` and hence loses its ability to print an error
message if assertions are deactivated.
"""
function require(packages::Vector{Symbol}; fun_name::String="",
                                           explanation::String="")
    check = any(isdefined(@__MODULE__, package) for package in packages)
    @assert check "no package from '$packages' loaded" *
        (fun_name == "" ? "" :
            " (one of them is required for executing `$fun_name`" *
            (explanation == "" ? "" : " " * explanation) * ")")
end
