import Base: *, ∈, isempty

export InverseLinearMap,
       an_element,
       constraints_list

struct InverseLinearMap{N<:Real, S<:LazySet{N},
                 NM, MAT<:AbstractMatrix{NM}} <: AbstractAffineMap{N, S}
    M::MAT
    X::S

    # default constructor with dimension match check
    function InverseLinearMap(M::MAT, X::S;
                                check_invertibility::Bool=true) where {N<:Real,
                                    S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}}
        @assert dim(X) == size(M, 1) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(X))"
        if check_invertibility
            @assert isinvertible(M) "the linear map is not invertible"
        end
        return new{N, S, NM, MAT}(M, X)
    end
end


# convenience constructor from a UniformScaling
function InverseLinearMap(M::UniformScaling{N}, X::LazySet) where {N}
    if M.λ == one(N)
        return X
    end
    return InverseLinearMap(Diagonal(fill(M.λ, dim(X))), X)
end


# convenience constructor from a scalar
function InverseLinearMap(α::Real, X::LazySet)
    n = dim(X)
    return InverseLinearMap(sparse(α * I, n, n), X)
end

# combine two linear maps into a single linear map
function InverseLinearMap(M::AbstractMatrix, lm::LinearMap)
    return InverseLinearMap(M * lm.M, lm.X)
end

# ZeroSet is "almost absorbing" for LinearMap (only the dimension changes)
function InverseLinearMap(M::AbstractMatrix{N}, Z::ZeroSet{N}) where {N}
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(Z))"
    return ZeroSet{N}(size(M, 1))
end

# EmptySet is absorbing for LinearMap
function InverseLinearMap(M::AbstractMatrix, ∅::EmptySet)
    return ∅
end


# --- AbstractAffineMap interface functions ---

function matrix(lm::InverseLinearMap)
    return lm.M
end

function vector(lm::InverseLinearMap{N}) where {N}
    return spzeros(N, dim(lm))
end

function set(lm::InverseLinearMap)
    return lm.X
end


function dim(lm::LinearMap)
    return size(lm.M, 1)
end

function ρ(d::AbstractVector, lm::InverseLinearMap)
    return ρ(transpose(lm.M) \ d, lm.X)
end

function σ(d::AbstractVector, lm::LinearMap)
    return lm.M \ σ(transpose(lm.M) \ d, lm.X)
end


```
Note that ``x ∈ M^{-1}⋅S`` iff ``M⋅x ∈ S``.
```
function ∈(x::AbstractVector, lm::InverseLinearMap)
    return lm.M * x ∈ lm.X
end



function an_element(lm::InverseLinearMap)
    return lm.M \ an_element(lm.X)
end




function vertices_list(lm::InverseLinearMap{N}; prune::Bool=true) where {N}
    # for a zero map, the result is just the list containing the origin
    if iszero(lm.M)
        return [zeros(N, dim(lm))]
    end

    # collect low-dimensional vertices lists
    vlist_X = vertices_list(lm.X)

    # create resulting vertices list
    vlist = Vector{Vector{N}}()
    sizehint!(vlist, length(vlist_X))
    for v in vlist_X
        push!(vlist, lm.M \ v)
    end

    return prune ? convex_hull(vlist) : vlist
end


function constraints_list(lm::InverseLinearMap{N}) where {N}
    return constraints_list(linear_map_inverse(lm.M, lm.X))
end


"""
    linear_map(M::AbstractMatrix{N}, lm::LinearMap{N}) where {N<:Real}

Return the linear map of a lazy linear map.

### Input

- `M`  -- matrix
- `lm` -- linear map

### Output

The polytope representing the linear map of the lazy linear map of a set.
"""
function linear_map(M::AbstractMatrix{N}, lm::InverseLinearMap{N}) where {N<:Real}
    return linear_map(M * inv(lm.M), lm.X) #TODO: agregar nota sobre calculo de inv(lm.M)
end

function concretize(lm::LinearMap)
    return linear_map(inv(lm.M), concretize(lm.X))
end
