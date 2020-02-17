eval(quote
    using .GeometryTypes: HyperRectangle, HyperCube, HyperSphere
end)

eval(load_geometry_types_conversions())
