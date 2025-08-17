const args1 = (1,)
const args12 = (1, 2)
const args123 = (1, 2, 3)
const args2 = (2,)
const global VALIDATE_DICT = Dict{Symbol,Tuple{Function,Any}}()

# unary set operations

# function validate_an_element(X::LazySet)
#     # require nonemptiness?
# end
# push!(VALIDATE_DICT, :an_element => (validate_an_element, args1))

function validate_area(X::LazySet)
    n = dim(X)
    if n != 2 && n != 3
        throw(DimensionMismatch("`area` requires a set of dimension 2 or 3 " *
                                "but received a $n-dimensional set"))
    end
    return true
end
push!(VALIDATE_DICT, :area => (validate_area, args1))

function validate_center(X::LazySet, i::Int)
    return validate_index(i, X; fun=center)
end
push!(VALIDATE_DICT, :center => (validate_center, args12))

# function validate_constraints_list(X::LazySet)
#     # require polyhedral set?
# end
# push!(VALIDATE_DICT, :constraints_list => (validate_constraints_list, args1))
# push!(VALIDATE_DICT, :constraints => (validate_constraints, args1))

function validate_diameter(p::Real)
    return validate_pnorm(p; fun=diameter)
end
push!(VALIDATE_DICT, :diameter => (validate_diameter, args2))

function validate_extrema(X::LazySet, i::Int)
    return validate_index(i, X; fun=extrema)
end
push!(VALIDATE_DICT, :extrema => (validate_extrema, args12))

function validate_high(X::LazySet, i::Int)
    return validate_index(i, X; fun=high)
end
push!(VALIDATE_DICT, :high => (validate_high, args12))

function validate_low(X::LazySet, i::Int)
    return validate_index(i, X; fun=low)
end
push!(VALIDATE_DICT, :low => (validate_low, args12))

function validate_norm(p::Real)
    return validate_pnorm(p; fun=norm)
end
push!(VALIDATE_DICT, :norm => (validate_norm, args2))

function validate_radius(p::Real)
    return validate_pnorm(p; fun=radius)
end
push!(VALIDATE_DICT, :radius => (validate_radius, args2))

function validate_radius_hyperrectangle(X::LazySet, i::Int)
    return validate_index(i, X; fun=radius_hyperrectangle)
end
push!(VALIDATE_DICT, :radius_hyperrectangle => (validate_radius_hyperrectangle, args12))

function validate_triangulate_faces(X::LazySet)
    n = dim(X)
    if n != 3
        throw(DimensionMismatch("`triangulate_faces` requires a set of dimension 3 " *
                                "but received a $n-dimensional set"))
    elseif !(ispolyhedral(X) && isbounded(X))
        throw(ArgumentError("`triangulate_faces` requires a polytopic set"))
    end
    return true
end
push!(VALIDATE_DICT, :triangulate_faces => (validate_triangulate_faces, args1))

# function validate_vertices_list(X::LazySet)
#     # require polytopic set?
# end
# push!(VALIDATE_DICT, :vertices_list => (validate_vertices_list, args1))
# push!(VALIDATE_DICT, :vertices => (validate_vertices, args1))

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

function validate_distance(x::AbstractVector, X::LazySet, p::Real)
    return validate_same_dim(x, X; fun=distance) && validate_pnorm(p; fun=distance)
end
push!(VALIDATE_DICT, :distance => (validate_distance, (1, 2, :p)))

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
    return validate_same_dim(x, X; fun=∈)
end
push!(VALIDATE_DICT, :∈ => (validate_in, args12))

function validate_is_interior_point(x::AbstractVector, X::LazySet, p::Real, ε::Number)
    if !validate_same_dim(x, X; fun=is_interior_point) || !validate_pnorm(p; fun=is_interior_point)
        return false
    elseif ε <= zero(ε)
        throw(ArgumentError("the tolerance must be strictly positive but is $ε"))
    end
    return true
end
push!(VALIDATE_DICT, :is_interior_point => (validate_is_interior_point, (1, 2, :p, :ε)))

function validate_linear_map(M::AbstractMatrix, X::LazySet)
    return validate_map_dim(M, X; fun=exponential_map)
end
push!(VALIDATE_DICT, :linear_map => (validate_linear_map, args12))

function validate_permute(X::LazySet, p::AbstractVector{Int})
    return validate_same_dim(p, X; fun=permute) && validate_index_vector(p, X; fun=permute)
end
push!(VALIDATE_DICT, :permute => (validate_permute, args12))

function validate_project(X::LazySet, block::AbstractVector{Int})
    n = dim(X)
    if isempty(block) || length(block) > n
        throw(DimensionMismatch("`project` requires an index vector for " *
                                "dimension $n but received a vector of length $(length(block))"))
    end
    return validate_index_vector(block, X; fun=project)
end
push!(VALIDATE_DICT, :project => (validate_project, args12))

# function validate_scale(α::Real, X::LazySet)
#     # require nonnegative factor?
# end
# push!(VALIDATE_DICT, :scale => (validate_scale, args12))
# push!(VALIDATE_DICT, :scale! => (validate_scale, args12))

function validate_support_function(d::AbstractVector, X::LazySet)
    return validate_same_dim(d, X; fun=ρ)
end
push!(VALIDATE_DICT, :ρ => (validate_support_function, args12))

function validate_support_vector(d::AbstractVector, X::LazySet)
    return validate_same_dim(d, X; fun=σ)
end
push!(VALIDATE_DICT, :σ => (validate_support_vector, args12))

function validate_translate(X::LazySet, v::AbstractVector)
    return validate_same_dim(v, X; fun=translate)
end
push!(VALIDATE_DICT, :translate => (validate_translate, args12))
push!(VALIDATE_DICT, :translate! => (validate_translate, args12))

# binary set operations

function validate_convex_hull(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=convex_hull)
end
push!(VALIDATE_DICT, :convex_hull => (validate_convex_hull, args12))

function validate_difference(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=difference)
end
push!(VALIDATE_DICT, :difference => (validate_difference, args12))

function validate_distance(X::LazySet, Y::LazySet, p::Real)
    return validate_same_dim(X, Y; fun=distance) && validate_pnorm(p; fun=distance)
end
# NOTE: dictionary entry was already added above for other `distance` method
# push!(VALIDATE_DICT, :distance => (validate_distance, (1, 2, :p)))

function validate_exact_sum(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=exact_sum)
end
push!(VALIDATE_DICT, :exact_sum => (validate_exact_sum, args12))

function validate_intersection(X::LazySet, Y::LazySet)
    return validate_same_dim(X, Y; fun=intersection)
end
push!(VALIDATE_DICT, :intersection => (validate_intersection, args12))
push!(VALIDATE_DICT, :intersection! => (validate_intersection, args12))

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
