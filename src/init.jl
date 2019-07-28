function __init__()
    @require Expokit = "a1e7a1ef-7a5d-5822-a38c-be74e1bb89f4" load_expokit()
    @require Optim = "429524aa-4258-5aef-a3af-852621145aeb" load_optim()
    @require Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029" load_polyhedra()
    @require Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f" load_distributions()
    @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" load_makie()
end

function load_expokit()
    eval(load_expokit_sparsematrixexp())
    eval(load_expokit_exponentialmap())
    eval(load_expokit_exponentialprojectionmap())
end

function load_optim()
    eval(load_optim_intersection())
end

function load_polyhedra()
    eval(load_polyhedra_abstractpolytope())
    eval(load_polyhedra_hpolytope())
    eval(load_polyhedra_hpolyhedron())
    eval(load_polyhedra_vpolytope())
end

function load_distributions()
    eval(load_distributions_samples())
end

function load_makie()
    eval(load_mesh())
end

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
