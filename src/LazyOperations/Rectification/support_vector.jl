"""
    σ(d::AbstractVector, R::Rectification)

Return a support vector of a rectification.

### Input

- `d` -- direction
- `R` -- rectification

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.
"""
@validate function σ(d::AbstractVector, R::Rectification)
    _compute_exact_representation!(R)
    return σ(d, R.cache.set)
end

"""
    σ(d::AbstractVector, R::Rectification{N, <:AbstractHyperrectangle}) where {N}

Return a support vector of the rectification of a hyperrectangular set.

### Input

- `d` -- direction
- `R` -- rectification of a hyperrectangular set

### Output

A support vector in the given direction.

### Algorithm

Let ``R(·)`` be the rectification of a vector respectively a set, and let ``H``
be a hyperrectangle. Then ``σ_{R(H)}(d) = R(σ_{H}(d))``.
"""
@validate function σ(d::AbstractVector, R::Rectification{N,<:AbstractHyperrectangle}) where {N}
    return rectify(σ(d, R.X))
end

"""
    σ(d::AbstractVector, R::Rectification{N, <:CartesianProduct}) where {N}

Return a support vector of the rectification of a Cartesian product of two sets.

### Input

- `d` -- direction
- `R` -- rectification of a Cartesian product of two sets

### Output

A support vector in the given direction.

### Notes

Note that this implementation creates new `Rectification` objects that do not
get preserved. Hence a second support-vector query does not benefit from the
computations in the first query. For this use case another implementation should
be added.

### Algorithm

Rectification distributes with the Cartesian product.
Let ``R(·)`` be the rectification of a set.
We can just query a support vector for ``R(X)`` and ``R(Y)`` recursively:
``σ_{R(X × Y)}(d) = σ_{R(X)}(d_X) × σ_{R(Y)}(d_Y)``, where ``x × y``
concatenates vectors ``x`` and ``y``.
"""
@validate function σ(d::AbstractVector, R::Rectification{N,<:CartesianProduct}) where {N}
    X, Y = first(R.X), second(R.X)
    n1 = dim(X)
    return vcat(σ(d[1:n1], Rectification(X)), σ(d[(n1 + 1):end], Rectification(Y)))
end

"""
    σ(d::AbstractVector,
      R::Rectification{N, <:CartesianProductArray}) where {N}

Return a support vector of the rectification of a Cartesian product of a
finite number of sets.

### Input

- `d` -- direction
- `R` -- rectification of a Cartesian product of a finite number of sets

### Output

A support vector in the given direction.

### Notes

Note that this implementation creates new `Rectification` objects that do not
get preserved. Hence a second support-vector query does not benefit from the
computations in the first query. For this use case another implementation should
be added.

### Algorithm

Rectification distributes with the Cartesian product.
Let ``R(·)`` be the rectification of a set.
We can just query a support vector for each subspace recursively:
``σ_{R(X_1 × ⋯ × X_m)}(d) = σ_{R(X_1)}(d_{X_1}) × ⋯ × σ_{R(X_m)}(d_{X_m})``,
where ``x × y`` concatenates vectors ``x`` and ``y``.
"""
@validate function σ(d::AbstractVector, R::Rectification{N,<:CartesianProductArray}) where {N}
    svec = similar(d)
    i = 1
    for X in R.X
        nX = dim(X)
        j = i + nX - 1
        svec[i:j] = σ(d[i:j], Rectification(X))
        i = j + 1
    end
    return svec
end
