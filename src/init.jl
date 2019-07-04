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
    eval(quote
    import .Polyhedra: polyhedron
    export polyhedron
    using .Polyhedra: HRep, VRep,
                      removehredundancy!, removevredundancy!
    import JuMP, GLPK

    function default_polyhedra_backend(P, N::Type{<:AbstractFloat})
        return Polyhedra.default_library(LazySets.dim(P), Float64)
    end

    function default_polyhedra_backend(P, N::Type{<:Rational})
        return Polyhedra.default_library(LazySets.dim(P), Rational{Int})
    end

    function default_lp_solver(N::Type{<:AbstractFloat})
        return JuMP.with_optimizer(GLPK.Optimizer)
    end

    function default_lp_solver(N::Type{<:Rational})
        return JuMP.with_optimizer(GLPK.Optimizer, method=:Exact)
    end
    end)

    eval(load_polyhedra_hpolytope())
    eval(load_polyhedra_hpolyhedron())
    eval(load_polyhedra_vpolytope())
    initialize_mesh()
end

function load_distributions()
    eval(load_distributions_samples())
end

function load_makie()
    initialize_mesh()
end

function initialize_mesh()
    if isdefined(@__MODULE__, :Polyhedra) &&
       isdefined(@__MODULE__, :Makie)

       eval(load_mesh())
    end
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
