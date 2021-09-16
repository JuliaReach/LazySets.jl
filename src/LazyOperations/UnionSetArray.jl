export UnionSetArray,
       array

# ========================================
# n-ary set union
# ========================================

"""
   UnionSetArray{N, S<:LazySet{N}}

Type that represents the set union of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

The union of convex sets is typically not convex.
"""
struct UnionSetArray{N, S<:LazySet{N}}
   array::Vector{S}
end

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

"""
   dim(cup::UnionSetArray)

Return the dimension of the set union of a finite number of sets.

### Input

- `cup` -- union of a finite number of sets

### Output

The ambient dimension of the union of a finite number of sets.
"""
function dim(cup::UnionSetArray)
   return dim(cup.array[1])
end

"""
   array(cup::UnionSetArray)

Return the array of the union of a finite number of sets.

### Input

- `cup` -- union of a finite number of sets

### Output

The array that holds the sets.
"""
function array(cup::UnionSetArray)
   return cup.array
end

"""
   σ(d::AbstractVector, cup::UnionSetArray; [algorithm]="support_vector")

Return the support vector of the union of a finite number of sets in a given
direction.

### Input

- `d`         -- direction
- `cup`       -- union of a finite number of sets
- `algorithm` -- (optional, default: "support_vector"): the algorithm to compute
                 the support vector; if "support_vector", use the support
                 vector of each argument; if "support_function" use the support
                 function of each argument and evaluate the support vector of
                 only one of them

### Output

The support vector in the given direction.

### Algorithm

The support vector of the union of a finite number of sets ``X₁, X₂, ...`` can
be obtained as the vector that maximizes the support function, i.e., it is
sufficient to find the ``\\argmax([ρ(d, X₂), ρ(d, X₂), ...])`` and evaluate its
support vector.

The default implementation, with option `algorithm="support_vector"`, computes
the support vector of all ``X₁, X₂, ...`` and then compares the support function
using a dot product.

If the support function can be computed more efficiently, the alternative
implementation `algorithm="support_function"` can be used, which evaluates the
support function of each set directly and then calls only the support vector of
one of the ``Xᵢ``.
"""
function σ(d::AbstractVector, cup::UnionSetArray; algorithm="support_vector")
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
   ρ(d::AbstractVector, cup::UnionSetArray)

Return the support function of the union of a finite number of sets in a given
direction.

### Input

- `d`   -- direction
- `cup` -- union of a finite number of sets

### Output

The support function in the given direction.

### Algorithm

The support function of the union of a finite number of sets ``X₁, X₂, ...``
can be obtained as the maximum of ``ρ(d, X₂), ρ(d, X₂), ...``.
"""
function ρ(d::AbstractVector, cup::UnionSetArray)
   A = array(cup)
   return maximum(Xi -> ρ(d, Xi), A)
end

"""
   an_element(cup::UnionSetArray)

Return some element of the union of a finite number of sets.

### Input

- `cup` -- union of a finite number of sets

### Output

An element in the union of a finite number of sets.

### Algorithm

We use `an_element` on the first wrapped set.
"""
function an_element(cup::UnionSetArray)
    for Xi in array(cup)
        if !isempty(Xi)
            return an_element(Xi)
        end
    end
    error("an empty set does not have any element")
end

"""
   ∈(x::AbstractVector, cup::UnionSetArray)

Check whether a given point is contained in the union of a finite number of
sets.

### Input

- `x`   -- point/vector
- `cup` -- union of a finite number of sets

### Output

`true` iff ``x ∈ cup``.
"""
function ∈(x::AbstractVector, cup::UnionSetArray)
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

Determine whether the union of a finite number of sets is bounded.

### Input

- `cup` -- union of a finite number of sets

### Output

`true` iff the union is bounded.
"""
function isbounded(cup::UnionSetArray)
   return all(isbounded, array(cup))
end

"""
    vertices_list(cup::UnionSetArray; apply_convex_hull::Bool=false,
                  backend=nothing)

Return the list of vertices of the union of a finite number of sets.

### Input

- `cup`               -- union of a finite number of sets
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
