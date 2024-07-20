function __init__()
    @require IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253" eval(load_IntervalArithmetic_convert())
    @require StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c" eval(load_genmat_hyperrectangle_static())
end
