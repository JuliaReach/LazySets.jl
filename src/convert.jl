import Base.convert

# convert methods for identity (no-ops)
for T in subtypes(LazySet, true)
    @eval begin
        Base.convert(::Type{$T}, X::$T) = X
    end
end

"""
    convert(::Type{VPolygon}, P::AbstractHPolygon)

Convert a polygon in constraint representation to a polygon in vertex
representation.

### Input

- `VPolygon` -- target type
- `P`        -- polygon in constraint representation

### Output

A polygon in vertex representation.
"""
function convert(::Type{VPolygon}, P::AbstractHPolygon)
    return tovrep(P)
end

# fast conversion from a 2D hyperrectangular set to a zonotope
function _convert_2D(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}
    c = center(H)
    rx = radius_hyperrectangle(H, 1)
    ry = radius_hyperrectangle(H, 2)
    G = _genmat_2D(c, rx, ry)
    return Zonotope(c, G)
end

@inline function _genmat_2D(c::AbstractVector{N}, rx, ry) where {N}
    flat_x = isapproxzero(rx)
    flat_y = isapproxzero(ry)
    ncols = !flat_x + !flat_y
    G = Matrix{N}(undef, 2, ncols)
    if !flat_x
        @inbounds begin
            G[1] = rx
            G[2] = zero(N)
        end
        if !flat_y
            @inbounds begin
                G[3] = zero(N)
                G[4] = ry
            end
        end
    elseif !flat_y
        @inbounds begin
            G[1] = zero(N)
            G[2] = ry
        end
    end
    return G
end

function load_genmat_2D_static()
    return quote
        @inline function _genmat_2D(c::SVector{L,N}, rx, ry) where {L,N}
            flat_x = isapproxzero(rx)
            flat_y = isapproxzero(ry)
            if !flat_x && !flat_y
                G = SMatrix{2,2,N,4}(rx, zero(N), zero(N), ry)
            elseif !flat_x && flat_y
                G = SMatrix{2,1,N,2}(rx, zero(N))
            elseif flat_x && !flat_y
                G = SMatrix{2,1,N,2}(zero(N), ry)
            else
                G = SMatrix{2,0,N,0}()
            end
            return G
        end

        # this function is type-stable but doesn't prune the generators according
        # to flat dimensions of H
        function _convert_2D_static(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}
            c = center(H)
            rx = radius_hyperrectangle(H, 1)
            ry = radius_hyperrectangle(H, 2)
            G = SMatrix{2,2,N,4}(rx, zero(N), zero(N), ry)
            return Zonotope(c, G)
        end

        function _convert_static(::Type{Zonotope},
                                 H::Hyperrectangle{N,<:SVector,<:SVector}) where {N}
            return Zonotope(center(H), _genmat_static(H))
        end
    end
end  # quote / load_genmat_2D_static

function convert(::Type{Zonotope}, H::AbstractHyperrectangle)
    dim(H) == 2 && return _convert_2D(Zonotope, H)
    return Zonotope(center(H), genmat(H))
end

function convert(::Type{Singleton},
                 cp::CartesianProduct{N,S1,S2}) where {N,S1<:AbstractSingleton,
                                                       S2<:AbstractSingleton}
    return Singleton(vcat(element(first(cp)), element(second(cp))))
end

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
    convert(::Type{Zonotope}, S::LinearMap{N, ZN})
        where {N, ZN<:AbstractZonotope}

Convert the lazy linear map of a zonotopic set to a zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- linear map of a zonotopic set

### Output

A zonotope.

### Algorithm

This method first applies the (concrete) linear map to the zonotopic set and
then converts the result to a `Zonotope` type.
"""
function convert(::Type{Zonotope}, S::LinearMap{N,ZN}) where {N,ZN<:AbstractZonotope}
    return convert(Zonotope, linear_map(S.M, S.X))
end

"""
    convert(::Type{Zonotope}, S::LinearMap{N, CartesianProduct{N, HN1, HN2}}
           ) where {N, HN1<:AbstractHyperrectangle,
                    HN2<:AbstractHyperrectangle}

Convert the lazy linear map of the Cartesian product of two hyperrectangular
sets to a zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- linear map of the Cartesian product of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method first converts the Cartesian product to a zonotope, and then
applies the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope},
                 S::LinearMap{N,CartesianProduct{N,HN1,HN2}}) where {N,HN1<:AbstractHyperrectangle,
                                                                     HN2<:AbstractHyperrectangle}
    return linear_map(S.M, convert(Zonotope, S.X))
end

"""
    convert(::Type{Zonotope},S::LinearMap{N, CartesianProductArray{N, HN}})
        where {N, HN<:AbstractHyperrectangle}

Convert the lazy linear map of the Cartesian product of a finite number of
hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- linear map of a `CartesianProductArray` of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method first converts the Cartesian product array to a zonotope, and then
applies the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope},
                 S::LinearMap{N,CartesianProductArray{N,HN}}) where {N,HN<:AbstractHyperrectangle}
    return linear_map(S.M, convert(Zonotope, S.X))
end

"""
    convert(::Type{Hyperrectangle}, cpa::CartesianProductArray{N, HN})
        where {N, HN<:AbstractHyperrectangle}

Convert the Cartesian product of a finite number of hyperrectangular sets to
a single hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `S`              -- Cartesian product array of hyperrectangular set

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `center` and `radius_hyperrectangle` methods of
`AbstractHyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 cpa::CartesianProductArray{N,HN}) where {N,HN<:AbstractHyperrectangle}
    n = dim(cpa)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    i = 1
    @inbounds for block_set in cpa
        j = i + dim(block_set) - 1
        c[i:j] = center(block_set)
        r[i:j] = radius_hyperrectangle(block_set)
        i = j + 1
    end
    return Hyperrectangle(c, r)
end

"""
    convert(::Type{Hyperrectangle}, cp::CartesianProduct{N, HN1, HN2})
        where {N, HN1<:AbstractHyperrectangle, HN2<:AbstractHyperrectangle}

Convert the Cartesian product of two hyperrectangular sets to a single
hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `S`              -- Cartesian product of two hyperrectangular sets

### Output

A hyperrectangle.

### Algorithm

The result is obtained by concatenating the center and radius of each
hyperrectangle. This implementation uses the `center` and
`radius_hyperrectangle` methods.
"""
function convert(::Type{Hyperrectangle},
                 cp::CartesianProduct{N,HN1,HN2}) where {N,HN1<:AbstractHyperrectangle,
                                                         HN2<:AbstractHyperrectangle}
    X, Y = first(cp), second(cp)
    c = vcat(center(X), center(Y))
    r = vcat(radius_hyperrectangle(X), radius_hyperrectangle(Y))
    return Hyperrectangle(c, r)
end

"""
    convert(::Type{Hyperrectangle},
            cpa::CartesianProductArray{N, IN}) where {N, IN<:Interval}

Convert the Cartesian product of a finite number of intervals to a single
hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `S`              -- Cartesian product array of intervals

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `min` and `max` methods of `Interval` to reduce
the allocations and improve performance (see LazySets#1143).
"""
function convert(::Type{Hyperrectangle},
                 cpa::CartesianProductArray{N,IN}) where {N,IN<:Interval}
    # since the sets are intervals, the dimension of cpa is its length
    n = length(array(cpa))
    l = Vector{N}(undef, n)
    h = Vector{N}(undef, n)
    @inbounds for (i, Ii) in enumerate(array(cpa))
        l[i] = min(Ii)
        h[i] = max(Ii)
    end
    return Hyperrectangle(; low=l, high=h)
end

"""
    convert(::Type{CartesianProduct{N, Interval{N}, Interval{N}}},
            H::AbstractHyperrectangle{N}) where {N}

Convert a two-dimensional hyperrectangle to the Cartesian product of two
intervals.

### Input

- `CartesianProduct` -- target type
- `H`                -- hyperrectangle

### Output

The Cartesian product of two intervals.
"""
function convert(::Type{CartesianProduct{N,Interval{N},Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N}
    @assert dim(H) == 2 "the hyperrectangle must be two-dimensional to " *
                        "convert it to the Cartesian product of two intervals, but it is " *
                        "$(dim(H))-dimensional; consider converting it to a " *
                        "`CartesianProductArray{$N, Interval{$N}}` instead"
    I1 = Interval(low(H, 1), high(H, 1))
    I2 = Interval(low(H, 2), high(H, 2))
    return CartesianProduct(I1, I2)
end

"""
    convert(::Type{CartesianProductArray{N, Interval{N}}},
            H::AbstractHyperrectangle{N}) where {N}

Convert a hyperrectangle to the Cartesian product array of intervals.

### Input

- `CartesianProductArray` -- target type
- `H`                     -- hyperrectangle

### Output

The Cartesian product of a finite number of intervals.
"""
function convert(::Type{CartesianProductArray{N,Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N}
    Iarray = [Interval(low(H, i), high(H, i)) for i in 1:dim(H)]
    return CartesianProductArray(Iarray)
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

"""
    convert(::Type{Hyperrectangle}, r::Rectification{N, AH})
        where {N, AH<:AbstractHyperrectangle}

Convert a rectification of a hyperrectangle to a hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `r`              -- rectification of a hyperrectangle

### Output

A `Hyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 r::Rectification{N,AH}) where {N,AH<:AbstractHyperrectangle}
    return rectify(r.X)
end

"""
    convert(::Type{Interval}, r::Rectification{N, IN}) where {N, IN<:Interval}

Convert a rectification of an interval to an interval.

### Input

- `Interval` -- target type
- `r`        -- rectification of an interval

### Output

An `Interval`.
"""
function convert(::Type{Interval},
                 r::Rectification{N,IN}) where {N,IN<:Interval}
    return Interval(rectify([min(r.X), max(r.X)]))
end

"""
    convert(::Type{MinkowskiSumArray},
            X::MinkowskiSum{N, ST, MinkowskiSumArray{N, ST}}) where {N, ST}

Convert the Minkowski sum of a Minkowski sum array to a Minkowski sum array.

### Input

- `MinkowskiSumArray`  -- target type
- `X`                  -- Minkowski sum of a Minkowski sum array

### Output

A Minkowski sum array.
"""
function convert(::Type{MinkowskiSumArray},
                 X::MinkowskiSum{N,ST,MinkowskiSumArray{N,ST}}) where {N,ST}
    return MinkowskiSumArray(vcat(first(X), array(second(X))))
end

"""
    convert(::Type{Interval}, ms::MinkowskiSum{N, IT, IT}) where {N, IT<:Interval}

Convert the Minkowski sum of two intervals to an interval.

### Input

- `Interval` -- target type
- `ms`       -- Minkowski sum of two intervals

### Output

An interval.
"""
function convert(::Type{Interval}, ms::MinkowskiSum{N,IT,IT}) where {N,IT<:Interval}
    return concretize(ms)
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

# make a copy of the constraints
convert(::Type{HPolytope}, P::HPolyhedron) = HPolytope(copy(constraints_list(P)))
convert(::Type{HPolyhedron}, P::HPolytope) = HPolyhedron(copy(constraints_list(P)))

for T in [HPolygon, HPolygonOpt, HPolytope, HPolyhedron]
    @eval begin
        function convert(::Type{$T}, P::Intersection)
            clist = vcat(constraints_list(first(P)), constraints_list(second(P)))
            return ($T)(clist)
        end

        function convert(::Type{$T}, P::IntersectionArray)
            clist = reduce(vcat, constraints_list.(array(P)))
            return ($T)(clist)
        end
    end
end

for T in subtypes(AbstractHPolygon, true)
    @eval begin
        """
            convert(::Type{$($T)}, X::LazySet; [check_boundedness]::Bool=true,
                    prune::Bool=true)

        Convert a two-dimensional polytopic set to a polygon in constraint
        representation.

        ### Input

        - `$($T)`             -- target type
        - `X`                 -- two-dimensional polytopic set
        - `check_boundedness` -- (optional, default `!isboundedtype(typeof(X))`) if
                                 `true` check whether the set `X` is bounded before
                                 creating the polygon
        - `prune`             -- (optional, default: `true`) flag for removing redundant
                                 constraints in the end

        ### Output

        A polygon in constraint representation.

        ### Algorithm

        We compute the list of constraints of `X`, then instantiate the polygon.
        """
        function convert(::Type{$T}, X::LazySet;
                         check_boundedness::Bool=!isboundedtype(typeof(X)),
                         prune::Bool=true)
            @assert dim(X) == 2 "set must be two-dimensional for conversion, but it " *
                                "has dimension $(dim(X))"
            if check_boundedness && !isbounded(X)
                throw(ArgumentError("expected a bounded set for conversion to `$T`"))
            end
            return $T(constraints_list(X); prune=prune)
        end

        """
            convert(T::Type{$($T)}, P::VPolygon)

        Convert a polygon in vertex representation to a polygon in constraint
        representation.

        ### Input

        - `$($T)` -- target type
        - `P`     -- polygon in vertex representation

        ### Output

        A polygon in constraint representation.
        """
        function convert(T::Type{$T}, P::VPolygon)
            return VPolygonModule.tohrep(P, T)
        end

        """
            convert(::Type{$($T)}, S::AbstractSingleton{N}) where {N}

        Convert a singleton to a polygon in constraint representation.

        ### Input

        - `$($T)` -- target type
        - `S`     -- singleton

        ### Output

        A polygon in constraint representation with the minimal number of constraints
        (three).
        """
        function convert(::Type{$T}, S::AbstractSingleton{N}) where {N}
            constraints_list = Vector{HalfSpace{N,Vector{N}}}(undef, 3)
            o = one(N)
            z = zero(N)
            v = element(S)
            constraints_list[1] = HalfSpace([o, o], v[1] + v[2])
            constraints_list[2] = HalfSpace([-o, z], -v[1])
            constraints_list[3] = HalfSpace([z, -o], -v[2])
            return $T(constraints_list)
        end

        """
            convert(::Type{$($T)}, L::LineSegment{N}) where {N}

        Convert a line segment to a polygon in constraint representation.

        ### Input

        - `$($T)` -- target type
        - `L`     -- line segment
        - `prune` -- (optional, default: `false`) flag for removing redundant
                     constraints in the end
        ### Output

        A flat polygon in constraint representation with the minimal number of
        constraints (four).
        """
        function convert(::Type{$T}, L::LineSegment{N}) where {N}
            H = $T{N}()
            c = halfspace_left(L.p, L.q)
            addconstraint!(H, c; prune=false)
            addconstraint!(H, HalfSpace(-c.a, -c.b); prune=false)
            line_dir = L.q - L.p
            c = HalfSpace(line_dir, dot(L.q, line_dir))
            addconstraint!(H, c; prune=false)
            line_dir = -line_dir
            addconstraint!(H, HalfSpace(line_dir, dot(L.p, line_dir)); prune=false)
            return H
        end
    end
end

"""
    convert(::Type{STAR}, P::AbstractPolyhedron{N}) where {N}

Convert a polyhedral set to a star set represented as a lazy affine map.

### Input

- `STAR` -- target type
- `P`    -- polyhedral set

### Output

A star set.
"""
function convert(::Type{STAR}, P::AbstractPolyhedron{N}) where {N}
    n = dim(P)
    c = zeros(N, n)
    V = Matrix(one(N) * I, n, n)
    return AffineMap(V, P, c)
end

"""
    convert(::Type{STAR}, X::Star)

Convert a star set to its equivalent representation as a lazy affine map.

### Input

- `STAR` -- target type
- `X`    -- star set

### Output

A star set.
"""
function convert(::Type{STAR}, X::Star)
    return AffineMap(X.V, X.P, X.c)
end

function convert(::Type{Hyperplane}, P::HPolyhedron; skip_check::Bool=false)
    # check that the number of constraints is fine
    if !skip_check && !ishyperplanar(P)
        throw(ArgumentError("the polyhedron is not hyperplanar: $P"))
    end

    # construct hyperplane from first constraint
    c1 = @inbounds first(constraints_list(P))
    return Hyperplane(c1.a, c1.b)
end

"""
    convert(::Type{SimpleSparsePolynomialZonotope}, Z::AbstractZonotope)

Convert a zonotope to a simple sparse polynomial zonotope.

### Input

- `SimpleSparsePolynomialZonotope` -- target type
- `Z`                              -- zonotopic set

### Output

A simple sparse polynomial zonotope.

### Algorithm

This method implements Proposition 3 in [1].

[1] Kochdumper, Althoff. *Sparse polynomial zonotopes - a novel set
representation for reachability analysis*. 2021
"""
function convert(::Type{SimpleSparsePolynomialZonotope}, Z::AbstractZonotope)
    c = center(Z)
    G = genmat(Z)
    n = ngens(Z)
    E = Matrix(1 * I, n, n)
    return SimpleSparsePolynomialZonotope(c, G, E)
end

"""
    convert(::Type{SimpleSparsePolynomialZonotope}, SPZ::SparsePolynomialZonotope)

Convert a sparse polynomial zonotope to simple sparse polynomial zonotope.

### Input

- `SimpleSparsePolynomialZonotope` -- target type
- `SPZ`                            -- sparse polynomial zonotope

### Output

A simple sparse polynomial zonotope.

### Algorithm

The method implements Proposition 3.1.4 from [1].

[1] Kochdumper, Niklas. *Extensions of polynomial zonotopes and their application to
verification of cyber-physical systems.* PhD diss., Technische Universität München, 2022.
"""
function convert(::Type{SimpleSparsePolynomialZonotope}, SPZ::SparsePolynomialZonotope)
    c = center(SPZ)
    G = hcat(genmat_dep(SPZ), genmat_indep(SPZ))
    n = ngens_indep(SPZ)
    E = cat(expmat(SPZ), Matrix(1 * I, n, n); dims=(1, 2))
    return SimpleSparsePolynomialZonotope(c, G, E)
end

"""
    convert(::Type{SparsePolynomialZonotope}, Z::AbstractZonotope{N}) where {N}

Convert a zonotope to sparse polynomial zonotope.

### Input

- `SparsePolynomialZonotope` -- target type
- `Z`                        -- zonotopic set

### Output

A sparse polynomial zonotope.

### Algorithm

The method implements Proposition 3.1.9 from [1].

[1] Kochdumper, Niklas. *Extensions of polynomial zonotopes and their application to
verification of cyber-physical systems.* PhD diss., Technische Universität München, 2022.
"""
function convert(::Type{SparsePolynomialZonotope}, Z::AbstractZonotope{N}) where {N}
    c = center(Z)
    G = genmat(Z)
    p = ngens(Z)
    E = Matrix(1 * I, p, p)
    GI = zeros(N, dim(Z), 0)
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

function convert(::Type{VPolytope}, T::Tetrahedron)
    return VPolytope(T.vertices)
end

function load_taylormodels_convert_polynomial_zonotope()
    return quote
        using .TaylorModels: TaylorModelN, polynomial, remainder, constant_term
        using .TaylorModels.TaylorSeries: coeff_table, set_variables, HomogeneousPolynomial,
                                          in_base, index_table, pos_table, TaylorN

        # implements Proposition 3.1.12 in thesis
        function convert(::Type{SparsePolynomialZonotope},
                         vTM::Vector{<:TaylorModelN{r,N}}) where {r,N}
            @assert _has_normalized_domain(vTM) "normalized domain (-1, 1) required"

            # upper bound on the number of terms/columns (ignores duplicates)
            # - 1 per iteration because we treat the constant terms separately
            num_coeffs = 0
            for TMi in vTM
                sum(TMi -> length(polynomial(TMi).coeffs) - 1, vTM)
            end

            n = length(vTM)
            c = Vector{N}(undef, n)
            Gs = Vector{Vector{N}}()
            GI_diag = Vector{N}(undef, n)
            Es = Vector{Vector{Int}}()

            total_columns = 0
            @inbounds for (i, TMi) in enumerate(vTM)
                pol = polynomial(TMi)
                rem = remainder(TMi)
                c[i] = IA.mid(rem) + constant_term(pol)
                GI_diag[i] = IA.radius(rem)
                for (order, term_order) in enumerate(pol.coeffs)
                    # we skip the first (= constant) term
                    if order == 1 || iszero(term_order)
                        continue
                    end
                    for (k, coeff_k) in enumerate(term_order.coeffs)
                        if iszero(coeff_k)
                            continue
                        end
                        Ej = coeff_table[order][k]
                        j = findfirst(e -> e == Ej, Es)
                        if isnothing(j)
                            total_columns += 1
                            j = total_columns
                            push!(Es, Ej)
                            push!(Gs, zeros(N, n))
                        end
                        Gs[j][i] += coeff_k
                    end
                end
            end
            G = stack(Gs)
            GI = remove_zero_columns(diagm(GI_diag))
            E = stack(Es)
            idx = uniqueID(n)
            return SparsePolynomialZonotope(c, G, GI, E, idx)
        end

        # implements Proposition 3.1.13 in thesis
        function convert(::Type{Vector{<:TaylorModelN}}, P::SparsePolynomialZonotope{N}) where {N}
            n = dim(P)
            p = nparams(P)  # number of parameters
            q = ngens_indep(P)  # independent/linear parameters
            h = ngens_dep(P)  # dependent/nonlinear parameters
            r = p + q
            c = center(P)
            G = genmat_dep(P)
            GI = genmat_indep(P)
            E = expmat(P)
            poly_order = polynomial_order(P)
            z = zeros(Int, q)
            # we need to rewrite the global variables
            set_variables("x"; order=poly_order, numvars=r)
            # the following vectors are shared for each polynomial
            rem = zero_itv(N)
            dom = sym_box(r, N)
            x0 = zero_box(r, N)

            vTM = Vector{TaylorModelN{r,N,N}}(undef, n)
            @inbounds for i in 1:n
                coeffs = zeros(HomogeneousPolynomial{N}, poly_order)

                # constant term
                coeffs[1] = c[i]

                # linear terms
                for j in 1:q
                    v = zeros(N, r)
                    v[p + j] = GI[i, j]
                    coeffs[2] += HomogeneousPolynomial(v, 1)
                end

                # nonlinear terms
                # TODO can be faster if loop is restructured so that index j is
                # created for each dimension
                for j in 1:h
                    Ej = E[:, j]
                    ord = sum(Ej)
                    G[i, j]
                    idx = pos_table[ord + 1][in_base(poly_order, Ej)]
                    l = length(coeff_table[ord + 1])
                    v = Vector(SingleEntryVector(idx, l, G[i, j]))
                    pj = HomogeneousPolynomial(v, ord)
                    coeffs[ord + 1] += pj
                end

                pol = TaylorN(coeffs, poly_order)
                vTM[i] = TaylorModelN(pol, rem, x0, dom)
            end
            return vTM
        end
    end
end  # quote / load_taylormodels_convert_polynomial_zonotope

function convert(::Type{Hyperrectangle}, Z::AbstractZonotope{N}) where {N}
    c = center(Z)
    n = length(c)
    r = zeros(N, n)
    @inbounds for cj in generators(Z)
        i = findfirst(!=(zero(N)), cj)
        if isnothing(i)
            continue
        end
        @assert isnothing(findfirst(!=(zero(N)), @view cj[(i + 1):end])) "the zonotope " *
                                                                         "is not hyperrectangular"
        r[i] += cj[i]  # `+` because to allow for multiple generators in dimension i
    end
    return Hyperrectangle(c, r)
end
