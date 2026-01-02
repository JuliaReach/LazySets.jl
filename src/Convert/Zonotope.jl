function load_StaticArraysCore_convert_Zonotope_Hyperrectangle_static()
    return quote
        function _convert_static(::Type{Zonotope},
                                 H::Hyperrectangle{N,<:SVector,<:SVector}) where {N}
            return Zonotope(center(H), _genmat_static(H))
        end
    end
end  # quote / load_StaticArraysCore_convert_Zonotope_Hyperrectangle_static

"""
    convert(::Type{Zonotope}, cp::CartesianProduct{N, HN1, HN2}) where {N,
            HN1<:AbstractHyperrectangle, HN2<:AbstractHyperrectangle}

Convert the Cartesian product of two hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- target type
- `cp`       -- Cartesian product of two hyperrectangular sets

### Output

This method falls back to the conversion of the Cartesian product to a single
hyperrectangle, and then from a hyperrectangle to a zonotope.
"""
function convert(::Type{Zonotope},
                 cp::CartesianProduct{N,HN1,HN2}) where {N,HN1<:AbstractHyperrectangle,
                                                         HN2<:AbstractHyperrectangle}
    return convert(Zonotope, convert(Hyperrectangle, cp))
end

"""
    convert(::Type{Zonotope}, cpa::CartesianProductArray{N, HN})
        where {N, HN<:AbstractHyperrectangle}

Convert the Cartesian product array of hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- target type
- `cpa`      -- Cartesian product array of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method falls back to the conversion of the Cartesian product to a single
hyperrectangle, and then from a hyperrectangle to a zonotope.
"""
function convert(::Type{Zonotope},
                 cpa::CartesianProductArray{N,HN}) where {N,HN<:AbstractHyperrectangle}
    return convert(Zonotope, convert(Hyperrectangle, cpa))
end

"""
    convert(::Type{Zonotope}, S::LinearMap{N, ZN}) where {N, ZN<:AbstractZonotope}

Convert the lazy linear map of a zonotopic set to a zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- linear map of a zonotopic set

### Output

A zonotope.

### Algorithm

This method first applies the (concrete) linear map to the zonotopic set and
then converts the result to a `Zonotope`.
"""
function convert(::Type{Zonotope}, S::LinearMap{N,ZN}) where {N,ZN<:AbstractZonotope}
    return convert(Zonotope, linear_map(S.M, S.X))
end

"""
    convert(::Type{Zonotope}, S::LinearMap{N, CartesianProduct{N, ZN1, ZN2}}
           ) where {N, ZN1<:AbstractZonotope, ZN2<:AbstractZonotope}

Convert the lazy linear map of the Cartesian product of two zonotopic
sets to a zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- linear map of the Cartesian product of zonotopic sets

### Output

A zonotope.

### Algorithm

This method first converts the Cartesian product to a zonotope and then applies
the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope},
                 S::LinearMap{N,CartesianProduct{N,ZN1,ZN2}}) where {N,ZN1<:AbstractZonotope,
                                                                     ZN2<:AbstractZonotope}
    return linear_map(S.M, convert(Zonotope, S.X))
end

"""
    convert(::Type{Zonotope},S::LinearMap{N, CartesianProductArray{N, ZN}})
        where {N, ZN<:AbstractZonotope}

Convert the lazy linear map of the Cartesian product of a finite number of
zonotopic sets to a zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- linear map of a `CartesianProductArray` of zonotopic sets

### Output

A zonotope.

### Algorithm

This method first converts the Cartesian product array to a zonotope and then
applies the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope},
                 S::LinearMap{N,CartesianProductArray{N,ZN}}) where {N,ZN<:AbstractZonotope}
    return linear_map(S.M, convert(Zonotope, S.X))
end

"""
    convert(::Type{Zonotope}, cp::CartesianProduct{N, ZN1, ZN2}
           ) where {N, ZN1<:AbstractZonotope, ZN2<:AbstractZonotope}

Convert the Cartesian product of two zonotopic sets to a new zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- Cartesian product of two zonotopic sets

### Output

A zonotope.

### Algorithm

The Cartesian product is obtained by:

- Concatenating the centers of each input zonotope.
- Arranging the generators in block-diagonal fashion, and filled with zeros in
  the off-diagonal; for this reason, the generator matrix of the returned
  zonotope is built as a sparse matrix.
"""
function convert(::Type{Zonotope},
                 cp::CartesianProduct{N,ZN1,ZN2}) where {N,ZN1<:AbstractZonotope,
                                                         ZN2<:AbstractZonotope}
    Z1, Z2 = first(cp), second(cp)
    c = vcat(center(Z1), center(Z2))
    G = blockdiag(sparse(genmat(Z1)), sparse(genmat(Z2)))
    return Zonotope(c, G)
end

function convert(::Type{Zonotope},
                 am::AbstractAffineMap{N,<:AbstractZonotope{N}}) where {N}
    Z1 = convert(Zonotope, linear_map(matrix(am), set(am)))
    translate!(Z1, vector(am))
    return Z1
end

"""
    convert(::Type{Zonotope}, cpa::CartesianProductArray{N, AZ})
        where {N, AZ<:AbstractZonotope}

Convert a Cartesian product array of zonotopic sets to a zonotope.

### Input

- `Zonotope` -- target type
- `cpa`       -- Cartesian product array of zonotopic sets

### Output

A zonotope with sparse matrix representation.
"""
function convert(::Type{Zonotope}, cpa::CartesianProductArray{N,AZ}) where {N,AZ<:AbstractZonotope}
    arr = array(cpa)
    c = reduce(vcat, center.(arr))
    G = reduce(blockdiag, sparse.(genmat.(arr)))
    return Zonotope(c, G)
end

function convert(::Type{Zonotope}, P::AbstractPolynomialZonotope)
    @assert iszero(ngens_dep(P)) "cannot convert a general polynomial zonotope to a Zonotope"

    return Zonotope(center(P), genmat_indep(P))
end
