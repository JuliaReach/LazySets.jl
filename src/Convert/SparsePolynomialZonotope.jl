"""
    convert(::Type{SparsePolynomialZonotope}, Z::AbstractZonotope; [algorithm]="GI")

Convert a zonotope to sparse polynomial zonotope.

### Input

- `SparsePolynomialZonotope` -- target type
- `Z`                        -- zonotopic set
- `algorithm`                -- (optional, default: `"GI"`) algorithm

### Output

A sparse polynomial zonotope.

### Algorithm

The `"GI"` method creates a polynomial zonotope with only independent generators.

The `"K21"` method creates a polynomial zonotope with only dependent generators,
implementing [Kochdumper21a; Proposition 3.1.9](@citet).
"""
function convert(::Type{SparsePolynomialZonotope}, Z::AbstractZonotope;
                 algorithm="GI")
    if algorithm == "GI"
        return _convert_SPZ_Z_GI(Z)
    elseif algorithm == "K21"
        return _convert_SPZ_Z_K21(Z)
    end
    throw(ArgumentError("invalid algorithm $algorithm"))
end

function _convert_SPZ_Z_K21(Z::AbstractZonotope{N}) where {N}
    c = center(Z)
    G = genmat(Z)
    p = ngens(Z)
    E = Matrix(1 * I, p, p)
    GI = zeros(N, dim(Z), 0)
    return SparsePolynomialZonotope(c, G, GI, E)
end

function _convert_SPZ_Z_GI(Z::AbstractZonotope{N}) where {N}
    c = center(Z)
    GI = genmat(Z)
    E = zeros(Int, 0, 0)
    G = zeros(N, length(c), 0)
    return SparsePolynomialZonotope(c, G, GI, E)
end

"""
    convert(::Type{SparsePolynomialZonotope}, SSPZ::SimpleSparsePolynomialZonotope)

Convert a simple sparse polynomial zonotope to a sparse polynomial zonotope.

### Input

- `SparsePolynomialZonotope` -- target type
- `SSPZ`                     -- simple sparse polynomial zonotope

### Output

A sparse polynomial zonotope.
"""
function convert(::Type{SparsePolynomialZonotope},
                 SSPZ::SimpleSparsePolynomialZonotope{N}) where {N}
    c = center(SSPZ)
    G = genmat(SSPZ)
    E = expmat(SSPZ)
    GI = Matrix{N}(undef, dim(SSPZ), 0)
    return SparsePolynomialZonotope(c, G, GI, E)
end
