function __init__()
    @require Expokit = "a1e7a1ef-7a5d-5822-a38c-be74e1bb89f4" include("init_Expokit.jl")
    @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("init_Makie.jl")
    @require Optim = "429524aa-4258-5aef-a3af-852621145aeb" include("init_Optim.jl")
    @require Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029" include("init_Polyhedra.jl")
    @require Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f" include("init_Distributions.jl")
    @require CDDLib = "3391f64e-dcde-5f30-b752-e11513730f60" eval(init_minkowski_sum_hrep())
end

function initialize_mesh()
    if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :Makie)
       eval(load_mesh())
    end
end

# dummy
function init_minkowski_sum_hrep() end

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
