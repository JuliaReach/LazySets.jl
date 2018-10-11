function __init__()
    @require Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029" load_polyhedra()
    @require Optim = "429524aa-4258-5aef-a3af-852621145aeb" load_optim()
end

function load_polyhedra()
    eval(load_polyhedra_abstractpolytope())
    eval(load_polyhedra_hpolytope())
    eval(load_polyhedra_hpolyhedron())
    eval(load_polyhedra_vpolytope())
    eval(load_polyhedra_concrete_intersection())
end

function load_optim()
    eval(load_optim_intersection())
end
