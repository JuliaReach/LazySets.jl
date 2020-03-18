function __init__()
    @require CDDLib = "3391f64e-dcde-5f30-b752-e11513730f60" include("Initialization/init_CDDLib.jl")
    @require Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f" include("Initialization/init_Distributions.jl")
    @require Expokit = "a1e7a1ef-7a5d-5822-a38c-be74e1bb89f4" include("Initialization/init_Expokit.jl")
    @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("Initialization/init_Makie.jl")
    @require MiniQhull = "1fd63fb3-3f00-53da-9a25-9f8c7140fe46" include("Initialization/init_MiniQhull.jl")
    @require Optim = "429524aa-4258-5aef-a3af-852621145aeb" include("Initialization/init_Optim.jl")
    @require Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029" include("Initialization/init_Polyhedra.jl")
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
