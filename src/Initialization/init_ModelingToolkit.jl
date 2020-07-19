eval(quote
    using .ModelingToolkit: get_variables,
                            gradient,
                            simplify,
                            substitute,
                            Operation

end)

eval(load_modeling_toolkit_hyperplane())
eval(load_modeling_toolkit_halfspace())
