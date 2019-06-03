function __init__()
    @require Expokit = "a1e7a1ef-7a5d-5822-a38c-be74e1bb89f4" load_expokit()
    @require Optim = "429524aa-4258-5aef-a3af-852621145aeb" load_optim()
    @require Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029" load_polyhedra()
    @require TaylorModels = "314ce334-5f6e-57ae-acf6-00b6e903104a" load_taylormodels()
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

function load_taylormodels()
    eval(load_taylormodels_overapproximation())
end
