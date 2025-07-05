const args1 = (1,)
const args12 = (1, 2)
const args123 = (1, 2, 3)
const global VALIDATE_DICT = Dict{Symbol,Tuple{Function,Tuple{Int,Vararg{Int}}}}()

# unary set operations

# function validate_an_element(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :an_element => (validate_an_element, args1))

function validate_area(X::LazySet)
    return validate_dims(X, (2, 3); fun=area)
end
push!(VALIDATE_DICT, :area => (validate_area, args1))

# function validate_center(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :center => (validate_center, args1))

# function validate_complement(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :complement => (validate_complement, args1))

# function validate_concretize(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :concretize => (validate_concretize, args1))

# function validate_constraints_list(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :constraints_list => (validate_constraints_list, args1))

# function validate_constraints(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :constraints => (validate_constraints, args1))

# function validate_convex_hull(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :convex_hull => (validate_convex_hull, args1))

# function validate_diameter(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :diameter => (validate_diameter, args1))

# function validate_dim(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :dim => (validate_dim, args1))

# function validate_eltype(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :eltype => (validate_eltype, args1))

# function validate_extrema(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :extrema => (validate_extrema, args1))

# function validate_high(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :high => (validate_high, args1))

# function validate_isbounded(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :isbounded => (validate_isbounded, args1))

# function validate_isempty(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :isempty => (validate_isempty, args1))

# function validate_isoperation(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :isoperation => (validate_isoperation, args1))

# function validate_ispolyhedral(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :ispolyhedral => (validate_ispolyhedral, args1))

# function validate_isuniversal(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :isuniversal => (validate_isuniversal, args1))

# function validate_low(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :low => (validate_low, args1))

# function validate_norm(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :norm => (validate_norm, args1))

# function validate_radius(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :radius => (validate_radius, args1))

# function validate_rectify(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :rectify => (validate_rectify, args1))

# function validate_vertices_list(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :vertices_list => (validate_vertices_list, args1))

# function validate_vertices(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :vertices => (validate_vertices, args1))

# function validate_volume(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :volume => (validate_volume, args1))

# mixed set operations

function validate_affine_map(M::AbstractMatrix, X::LazySet, v::AbstractVector)
    m = size(M, 1)
    if m != length(v)
        throw(DimensionMismatch("`affine_map` requires compatible matrix and vector dimensions " *
                                "but received a matrix of size $m and a vector of size " *
                                "$(length(v))"))
    end
    return validate_map_dim(M, X; fun=affine_map)
end
push!(VALIDATE_DICT, :affine_map => (validate_affine_map, args123))

function validate_distance(x::AbstractVector, X::LazySet)
    return validate_same_dims(x, X; fun=distance)
end
push!(VALIDATE_DICT, :distance => (validate_distance, args12))

function validate_exponential_map(M::AbstractMatrix, X::LazySet)
    m, n = size(M)
    if m != n
        throw(DimensionMismatch("`exponential_map` requires a quadratic matrix " *
                                "but received a matrix of size $m × $n"))
    end
    return validate_map_dim(M, X; fun=exponential_map)
end
push!(VALIDATE_DICT, :exponential_map => (validate_exponential_map, args12))

function validate_in(x::AbstractVector, X::LazySet)
    return validate_same_dims(x, X; fun=∈)
end
push!(VALIDATE_DICT, :∈ => (validate_in, args12))

function validate_is_interior_point(x::AbstractVector, X::LazySet)
    return validate_same_dims(x, X; fun=is_interior_point)
end
push!(VALIDATE_DICT, :is_interior_point => (validate_is_interior_point, args12))

function validate_linear_map(M::AbstractMatrix, X::LazySet)
    return validate_map_dim(M, X; fun=exponential_map)
end
push!(VALIDATE_DICT, :linear_map => (validate_linear_map, args12))

function validate_permute(X::LazySet, p::AbstractVector{Int})
    return validate_same_dims(p, X; fun=permute) && validate_index_vector(p, X)
end
push!(VALIDATE_DICT, :permute => (validate_permute, args12))

function validate_project(X::LazySet, block::AbstractVector{Int})
    return validate_index_vector_length(block, X) && validate_index_vector(block, X; fun=project)
end
push!(VALIDATE_DICT, :project => (validate_project, args12))

# function validate_sample(X::LazySet)
#     # nothing to validate
# end
# push!(VALIDATE_DICT, :sample => (validate_sample, args1))

function validate_scale(α::AbstractVector, X::LazySet)
    return validate_same_dims(α, X; fun=scale)
end
push!(VALIDATE_DICT, :scale => (validate_scale, args12))
push!(VALIDATE_DICT, :scale! => (validate_scale, args12))

function validate_support_function(d::AbstractVector, X::LazySet)
    return validate_same_dims(d, X; fun=ρ)
end
push!(VALIDATE_DICT, :ρ => (validate_support_function, args12))

function validate_support_vector(d::AbstractVector, X::LazySet)
    return validate_same_dims(d, X; fun=σ)
end
push!(VALIDATE_DICT, :σ => (validate_support_vector, args12))

function validate_translate(X::LazySet, v::AbstractVector)
    return validate_same_dims(v, X; fun=translate)
end
push!(VALIDATE_DICT, :translate => (validate_translate, args12))
push!(VALIDATE_DICT, :translate! => (validate_translate, args12))

# binary set operations

function validate_cartesian_product(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=cartesian_product)
end
push!(VALIDATE_DICT, :cartesian_product => (validate_cartesian_product, args12))

function validate_convex_hull(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=convex_hull)
end
push!(VALIDATE_DICT, :convex_hull => (validate_convex_hull, args12))

function validate_difference(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=difference)
end
push!(VALIDATE_DICT, :difference => (validate_difference, args12))

function validate_distance(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=distance)
end
push!(VALIDATE_DICT, :distance => (validate_distance, args12))

function validate_exact_sum(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=exact_sum)
end
push!(VALIDATE_DICT, :exact_sum => (validate_exact_sum, args12))

function validate_intersection(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=intersection)
end
push!(VALIDATE_DICT, :intersection => (validate_intersection, args12))

function validate_isdisjoint(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=isdisjoint)
end
push!(VALIDATE_DICT, :isdisjoint => (validate_isdisjoint, args12))

function validate_isequivalent(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=isequivalent)
end
push!(VALIDATE_DICT, :isequivalent => (validate_isequivalent, args12))

function validate_isstrictsubset(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=⊂)
end
push!(VALIDATE_DICT, :⊂ => (validate_isstrictsubset, args12))

function validate_issubset(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=⊆)
end
push!(VALIDATE_DICT, :⊆ => (validate_issubset, args12))

function validate_linear_combination(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=linear_combination)
end
push!(VALIDATE_DICT, :linear_combination => (validate_linear_combination, args12))

function validate_minkowski_difference(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=minkowski_difference)
end
push!(VALIDATE_DICT, :minkowski_difference => (validate_minkowski_difference, args12))

function validate_minkowski_sum(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=minkowski_sum)
end
push!(VALIDATE_DICT, :minkowski_sum => (validate_minkowski_sum, args12))
