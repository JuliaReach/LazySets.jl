export UnionSetArray,
       array

"""
   UnionSetArray{N, S<:LazySet{N}} <: LazySet{N}

Type that represents the set union of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

The union of convex sets is typically not convex.
"""
struct UnionSetArray{N,S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

"""
    UnionSet!(X, Y)

Convenience function to compute the lazy union and modify `UnionSetArray`s
in-place.
"""
function UnionSet! end

isoperationtype(::Type{<:UnionSetArray}) = true
isconvextype(::Type{<:UnionSetArray}) = false

# constructor for an empty union with optional size hint and numeric type
function UnionSetArray(n::Int=0, N::Type=Float64)
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return UnionSetArray(arr)
end

# EmptySet is the neutral element for UnionSetArray
@neutral(UnionSetArray, EmptySet)

# Universe is the absorbing element for UnionSetArray
@absorbing(UnionSetArray, Universe)

# add functions connecting UnionSet and UnionSetArray
@declare_array_version(UnionSet, UnionSetArray)

concretize(cup::UnionSetArray) = UnionSetArray([concretize(X) for X in array(cup)])

"""
   dim(cup::UnionSetArray)

Return the dimension of the union of a finite number of sets.

### Input

- `cup` -- union of a finite number of sets

### Output

The ambient dimension of the union of a finite number of sets, or `0` if there
is no set in the array.
"""
function dim(cup::UnionSetArray)
    return length(cup.array) == 0 ? 0 : dim(cup.array[1])
end

"""
   array(cup::UnionSetArray)

Return the array of the union of a finite number of sets.

### Input

- `cup` -- union of a finite number of sets

### Output

The array of the union.
"""
function array(cup::UnionSetArray)
    return cup.array
end

"""
   σ(d::AbstractVector, cup::UnionSetArray; [algorithm]="support_vector")

Return a support vector of the union of a finite number of sets in a given
direction.

### Input

- `d`         -- direction
- `cup`       -- union of a finite number of sets
- `algorithm` -- (optional, default: "support_vector"): the algorithm to compute
                 the support vector; if "support_vector", use the support
                 vector of each argument; if "support_function", use the support
                 function of each argument and evaluate the support vector of
                 only one of them

### Output

A support vector in the given direction.

### Algorithm

The support vector of the union of a finite number of sets ``X₁, X₂, ...`` can
be obtained as the vector that maximizes the support function, i.e., it is
sufficient to find the ``\\argmax([ρ(d, X₂), ρ(d, X₂), ...])`` and evaluate its
support vector.

The default implementation, with option `algorithm="support_vector"`, computes
the support vector of all ``X₁, X₂, ...`` and then compares the support function
using the dot product.

If the support function can be computed more efficiently, the alternative
implementation `algorithm="support_function"` can be used, which evaluates the
support function of each set directly and then calls only the support vector of
one of the ``Xᵢ``.
"""
@validate function σ(d::AbstractVector, cup::UnionSetArray; algorithm="support_vector")
    arr = array(cup)

    if algorithm == "support_vector"
        return _σ_union(d, arr)

    elseif algorithm == "support_function"
        m = argmax(i -> ρ(d, @inbounds arr[i]), eachindex(arr))
        return σ(d, arr[m])

    else
        throw(ArgumentError("algorithm $algorithm unknown"))
    end
end

function _σ_union(d::AbstractVector, sets)
    σmax = d
    N = eltype(d)
    ρmax = N(-Inf)
    for Xi in sets
        σX = σ(d, Xi)
        ρX = dot(d, σX)
        if ρX > ρmax
            ρmax = ρX
            σmax = σX
        end
    end
    return σmax
end

"""
   ρ(d::AbstractVector, cup::UnionSetArray)

Evaluate the support function of the union of a finite number of sets in a given
direction.

### Input

- `d`   -- direction
- `cup` -- union of a finite number of sets

### Output

The evaluation of the support function in the given direction.

### Algorithm

The support function of the union of a finite number of sets ``X₁, X₂, ...``
can be obtained as the maximum of ``ρ(d, X₂), ρ(d, X₂), ...``.
"""
@validate function ρ(d::AbstractVector, cup::UnionSetArray)
    return maximum(Xi -> ρ(d, Xi), array(cup))
end

"""
   an_element(cup::UnionSetArray)

Return some element of the union of a finite number of sets.

### Input

- `cup` -- union of a finite number of sets

### Output

An element in the union of a finite number of sets.

### Algorithm

We use `an_element` on the first non-empty wrapped set.
"""
function an_element(cup::UnionSetArray)
    for Xi in cup
        if !isempty(Xi)
            return an_element(Xi)
        end
    end
    return throw(ArgumentError("an empty set does not have any element"))
end

"""
    in(x::AbstractVector, cup::UnionSetArray)

Check whether a given point is contained in the union of a finite number of
sets.

### Input

- `x`   -- point/vector
- `cup` -- union of a finite number of sets

### Output

`true` iff ``x ∈ cup``.
"""
@validate function in(x::AbstractVector, cup::UnionSetArray)
    return any(X -> x ∈ X, array(cup))
end

"""
   isempty(cup::UnionSetArray)

Check whether the union of a finite number of sets is empty.

### Input

- `cup` -- union of a finite number of sets

### Output

`true` iff the union is empty.
"""
function isempty(cup::UnionSetArray)
    return all(isempty, array(cup))
end

"""
   isbounded(cup::UnionSetArray)

Check whether the union of a finite number of sets is bounded.

### Input

- `cup` -- union of a finite number of sets

### Output

`true` iff the union is bounded.
"""
function isbounded(cup::UnionSetArray)
    return all(isbounded, array(cup))
end

function isboundedtype(::Type{<:UnionSetArray{N,S}}) where {N,S}
    return isboundedtype(S)
end

"""
    vertices_list(cup::UnionSetArray; [apply_convex_hull]::Bool=false,
                  [backend]=nothing)

Return a list of vertices of the union of a finite number of sets.

### Input

- `cup`               -- union of a finite number of sets
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)

### Output

A list of vertices, possibly reduced to the list of vertices of the convex hull.
"""
function vertices_list(cup::UnionSetArray;
                       apply_convex_hull::Bool=false,
                       backend=nothing)
    vlist = vcat([vertices_list(Xi) for Xi in cup]...)
    if apply_convex_hull
        convex_hull!(vlist; backend=backend)
    end
    return vlist
end

@validate function linear_map(M::AbstractMatrix, cup::UnionSetArray)
    return UnionSetArray([linear_map(M, X) for X in cup])
end

@validate function project(cup::UnionSetArray, block::AbstractVector{Int}; kwargs...)
    return UnionSetArray([project(X, block; kwargs...) for X in cup])
end

@validate function translate(cup::UnionSetArray, v::AbstractVector)
    return UnionSetArray([translate(X, v) for X in cup])
end
