export UnionSetArray,
       array

# ========================================
# n-ary set union
# ========================================

"""
   UnionSetArray{N<:Real, S<:LazySet{N}}

Type that represents the set union of a finite number of convex sets.

### Fields

- `array` -- array of convex sets
"""
struct UnionSetArray{N<:Real, S<:LazySet{N}}
   array::Vector{S}
end

isoperationtype(::Type{<:UnionSetArray}) = true
isconvextype(::Type{<:UnionSetArray}) = false

# EmptySet is the neutral element for UnionSetArray
@neutral(UnionSetArray, EmptySet)

# Universe is the absorbing element for UnionSetArray
@absorbing(UnionSetArray, Universe)

# add functions connecting UnionSet and UnionSetArray
@declare_array_version(UnionSet, UnionSetArray)

"""
   dim(cup::UnionSetArray)

Return the dimension of the set union of a finite number of convex sets.

### Input

- `cup` -- union of a finite number of convex sets

### Output

The ambient dimension of the union of a finite number of convex sets.
"""
function dim(cup::UnionSetArray)
   return dim(cup.array[1])
end

"""
   array(cup::UnionSetArray{N, S}) where {N<:Real, S<:LazySet{N}}

Return the array of a union of a finite number of convex sets.

### Input

- `cup` -- union of a finite number of convex sets

### Output

The array that holds the union of a finite number of convex sets.
"""
function array(cup::UnionSetArray{N, S}) where {N<:Real, S<:LazySet{N}}
   return cup.array
end

"""
   σ(d::AbstractVector{N}, cup::UnionSetArray{N}; [algorithm]="support_vector") where {N<:Real}

Return the support vector of the union of a finite number of convex sets in
a given direction.

### Input

- `d`         -- direction
- `cup`       -- union of a finite number of convex sets
- `algorithm` -- (optional, default: "support_vector"): the algorithm to compute
                the support vector; if "support_vector", use the support
                vector of each argument; if "support_function" use the support
                function of each argument and evaluate the support vector of only
                one of them

### Output

The support vector in the given direction.

### Algorithm

The support vector of the union of a finite number of convex sets ``X₁, X₂, ...``
can be obtained as the vector that maximizes the support function, i.e.
it is sufficient to find the ``\\argmax(ρ(d, X₂), ρ(d, X₂), ...])`` and evaluate
its support vector.

The default implementation, with option `algorithm="support_vector"`, computes
the support vector of all ``X₁, X₂, ...`` and then compares the support function using
a dot product. If it happens that the support function can be more efficiently
computed (without passing through the support vector), consider using the alternative
`algorithm="support_function"` implementation, which evaluates the support function
of each set directly and then calls only the support vector of one of the ``Xᵢ``.
"""
function σ(d::AbstractVector{N}, cup::UnionSetArray{N}; algorithm="support_vector") where {N<:Real}
   A = array(cup)
   if algorithm == "support_vector"
       σarray = map(Xi -> σ(d, Xi), A)
       ρarray = map(vi -> dot(d, vi), σarray)
       m = argmax(ρarray)
       return σarray[m]
   elseif algorithm == "support_function"
       ρarray = map(Xi -> ρ(d, Xi), A)
       m = argmax(ρarray)
       return σ(d, A[m])
   else
       error("algorithm $algorithm for the support vector of a `UnionSetArray` is unknown")
   end
end

"""
   ρ(d::AbstractVector{N}, cup::UnionSetArray{N}) where {N<:Real}

Return the support function of the union of a finite number of convex sets in
a given direction.

### Input

- `d`   -- direction
- `cup` -- union of a finite number of convex sets

### Output

The support function in the given direction.

### Algorithm

The support function of the union of a finite number of convex sets ``X₁, X₂, ...``
can be obtained as the maximum of ``ρ(d, X₂), ρ(d, X₂), ...``.
"""
function ρ(d::AbstractVector{N}, cup::UnionSetArray{N}) where {N<:Real}
   A = array(cup)
   ρarray = map(Xi -> ρ(d, Xi), A)
   return maximum(ρarray)
end

"""
   an_element(cup::UnionSetArray{N}) where {N<:Real}

Return some element of a union of a finite number of convex sets.

### Input

- `cup` -- union of a finite number of convex sets

### Output

An element in the union of a finite number of convex sets.

### Algorithm

We use `an_element` on the first wrapped set.
"""
function an_element(cup::UnionSetArray{N}) where {N<:Real}
   return an_element(array(cup)[1])
end

"""
   ∈(x::AbstractVector{N}, cup::UnionSetArray{N}) where {N<:Real}

Check whether a given point is contained in a union of a finite number of convex
sets.

### Input

- `x`   -- point/vector
- `cup` -- union of a finite number of convex sets

### Output

`true` iff ``x ∈ cup``.
"""
function ∈(x::AbstractVector{N}, cup::UnionSetArray{N}) where {N<:Real}
   return any(X -> x ∈ X, array(cup))
end

"""
   isempty(cup::UnionSetArray)

Check whether a union of a finite number of convex sets is empty.

### Input

- `cup` -- union of a finite number of convex sets

### Output

`true` iff the union is empty.
"""
function isempty(cup::UnionSetArray)
   return all(isempty, array(cup))
end

"""
   isbounded(cup::UnionSetArray)

Determine whether a union of a finite number of convex sets is bounded.

### Input

- `cup` -- union of a finite number of convex sets

### Output

`true` iff the union is bounded.
"""
function isbounded(cup::UnionSetArray)
   return all(isbounded, array(cup))
end

"""
    vertices_list(cup::UnionSetArray; apply_convex_hull::Bool=false,
                  backend=nothing)

Return the list of vertices of a union of a finite number of convex sets.

### Input

- `cup`               -- union of a finite number of convex sets
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)

### Output

The list of vertices, possibly reduced to the list of vertices of the convex
hull.
"""
function vertices_list(cup::UnionSetArray;
                       apply_convex_hull::Bool=false,
                       backend=nothing)
    vlist = vcat([vertices_list(Xi) for Xi in array(cup)]...)
    if apply_convex_hull
        convex_hull!(vlist, backend=backend)
    end
    return vlist
end
