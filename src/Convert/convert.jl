import Base.convert

# convert methods for identity (no-ops)
for T in subtypes(LazySet, true)
    @eval begin
        Base.convert(::Type{$T}, X::$T) = X
    end
end

# =======================================
# `convert` methods sorted by target type
# =======================================

include("AffineMap.jl")
include("CartesianProduct.jl")
include("CartesianProductArray.jl")
include("DensePolynomialZonotope.jl")
include("HPolyhedron.jl")
include("HPolytope.jl")
include("Hyperplane.jl")
include("Hyperrectangle.jl")
include("Interval.jl")
include("MinkowskiSumArray.jl")
include("SimpleSparsePolynomialZonotope.jl")
include("Singleton.jl")
include("SparsePolynomialZonotope.jl")
include("TaylorModelN.jl")
include("VPolygon.jl")
include("VPolytope.jl")
include("Zonotope.jl")

# ===============================================
# `convert` methods with multiple by target types
# ===============================================

for T in [HPolygon, HPolygonOpt, HPolytope, HPolyhedron]
    @eval begin
        function Base.convert(::Type{$T}, P::Intersection)
            clist = vcat(constraints_list(first(P)), constraints_list(second(P)))
            return ($T)(clist)
        end

        function Base.convert(::Type{$T}, P::IntersectionArray)
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
        function Base.convert(::Type{$T}, X::LazySet;
                              check_boundedness::Bool=!isboundedtype(typeof(X)),
                              prune::Bool=true)
            @assert dim(X) == 2 "set must be two-dimensional for conversion, but it " *
                                "has dimension $(dim(X))"
            if check_boundedness && !isbounded(X)
                throw(ArgumentError("expected a bounded set for conversion to `$($T)`"))
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
        function Base.convert(T::Type{$T}, P::VPolygon)
            return VPolygonModule.tohrep(P, T)
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
