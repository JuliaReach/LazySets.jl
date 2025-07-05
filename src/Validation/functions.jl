const args1 = (1,)
const args12 = (1, 2)
const global VALIDATE_DICT = Dict{Symbol,Tuple{Function,Tuple{Int}}}()

# unary set operations

function validate_area(X::LazySet)
    return validate_dims(X, (2, 3); fun=area)
end
push!(VALIDATE_DICT, :area => (validate_area, args1))

# TODO add unary functions

# mixed set operations

function validate_affine_map()

end
# push!(VALIDATE_DICT, : => (validate_, args1))

function validate_distance(x::AbstractVector, X::LazySet)

end

function validate_exponential_map()

end

function validate_in()

end

function validate_is_interior_point()

end

function validate_linear_map()

end

function validate_permute()

end

function validate_project()

end

function validate_sample()

end

function validate_scale()

end

function validate_support_function()

end

function validate_support_vector()

end

function validate_translate()

end

# binary set operations

function validate_convex_hull(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=convex_hull)
end

function validate_cartesian_product()

end

function validate_difference()

end

function validate_distance(X::LazySet, Y::LazySet)

end

function validate_exact_sum()

end

function validate_intersection()

end

function validate_isdisjoint()

end

function validate_isequivalent()

end

function validate_isstrictsubset()

end

function validate_issubset()

end

function validate_linear_combination()

end

function validate_minkowski_difference()

end

function validate_minkowski_sum()

end
