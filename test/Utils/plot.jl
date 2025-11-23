using LazySets, Test, LinearAlgebra, SparseArrays
import LazySets.RecipesBase
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

# define `plot` function as `RecipesBase.apply_recipe`
struct DummyBackend <: RecipesBase.AbstractBackend end
struct DummyPlot <: RecipesBase.AbstractPlot{DummyBackend} end
Base.length(::DummyPlot) = 0
dict = Dict{Symbol,Any}(:plot_object => DummyPlot())
plot(args...; kwargs...) = RecipesBase.apply_recipe(dict, args...; kwargs...)
RecipesBase.is_key_supported(s::Symbol) = true

for N in @tN([Float64, Float32, Rational{Int}])
    p0 = zero(N)
    p1 = one(N)
    for n in [1, 2]
        v0 = zeros(N, n)
        v1 = ones(N, n)

        # ---------------
        # set definitions
        # ---------------

        # bounded basic set types
        es = EmptySet{N}(n)
        b1 = Ball1(v0, p1)
        bi = BallInf(v0, p1)
        hr = Hyperrectangle(v0, v1)
        st = Singleton(v1)
        zt = Zonotope(v0, Diagonal(ones(N, n)))
        ztmd = ZonotopeMD(v0, zeros(N, n, n), ones(N, n))
        zs = ZeroSet{N}(n)
        itv = Interval(p0, p1)
        if n == 2
            ls = LineSegment(v0, v1)
        end
        if N <: AbstractFloat
            hpa = convert(HParallelotope, hr)
        else
            @test_broken convert(HParallelotope, hr) isa HParallelotope  # TODO `convert` not available for non-float
            hpa = HParallelotope(n == 1 ? hcat(N[1]) : N[0 -1; 1 0], ones(N, 2 * n))
        end

        # -------------------------------------------
        # plot polytopes
        # -------------------------------------------
        plot(es)
        plot(b1)
        plot(bi)
        plot(hr)
        plot(st)
        plot(zt)
        if n == 2
            plot(ztmd)  # TODO crashes in 1D
        end
        plot(zs)
        if n == 1
            plot(itv)
        elseif n == 2
            plot(ls)
        end
        if N <: AbstractFloat
            plot(hpa)
        end

        if N == Rational{Int}
            # rationals do not support epsilon-close approximation
            continue
        end

        # bounded basic set types which are not polytopes
        b2 = Ball2(v0, p1)
        bp = Ballp(N(1.5), v0, p1)
        el = Ellipsoid(v0, Diagonal(ones(N, n)))
        dpz = convert(DensePolynomialZonotope, zt)
        spz = convert(SparsePolynomialZonotope, zt)
        sspz = convert(SimpleSparsePolynomialZonotope, zt)
        if n == 2
            ncp = Polygon([v0, v1])
        end

        # unary set operations
        spI = SparseMatrixCSC{N}(2I, n, n)
        smIe = SparseMatrixExp(spI)
        em = ExponentialMap(smIe, b2)
        psme = ProjectionSparseMatrixExp(spI, smIe, spI)
        epm = ExponentialProjectionMap(psme, b2)

        # binary set operations
        its = Intersection(b1, bi)
        itsa = IntersectionArray([b1, bi])

        # polygon/polytope types
        constraints = [LinearConstraint([p1, p0], p1),
                       LinearConstraint([p0, p1], p1),
                       LinearConstraint([-p1, p0], p0),
                       LinearConstraint([p0, -p1], p0)]
        hpg = HPolygon(constraints)
        hpgo = HPolygonOpt(constraints)
        hpt = HPolytope(constraints)
        hpt_empty = HPolytope([HalfSpace(N[1, 0], N(0)), HalfSpace(N[-1, 0], N(-1)),
                               HalfSpace(N[0, 1], N(0)), HalfSpace(N[0, -1], N(0))])
        vlist = vertices_list(bi)
        if n == 2
            vpg = VPolygon(vlist)
        end
        vpt = VPolytope(vlist)

        # unbounded sets
        hs = HalfSpace(v1, p1)
        hp = Hyperplane(v1, p1)
        hph = HPolyhedron(constraints)
        uni = Universe{N}(n)
        l = Line(v0, v1)
        if n == 2
            l2D = Line2D(v1, p1)
        end
        sta = convert(Star, bi)

        # unary set operations
        sih = SymmetricIntervalHull(b1)
        if n == 1
            lm = LinearMap(hcat(N[2]), bi)
        elseif n == 2
            lm = LinearMap(N[2 1; 1 2], bi)
        end
        rm = ResetMap(bi, Dict(1 => N(1)))

        # binary set operations
        ch = ConvexHull(b1, bi)
        cha = ConvexHullArray([b1, bi])
        ms = MinkowskiSum(b1, bi)
        msa = MinkowskiSumArray([b1, bi])
        cms = CachedMinkowskiSumArray([b1, bi])
        if n == 2
            cp = CartesianProduct(itv, itv)
            cpa = CartesianProductArray([itv, itv])
        end
        us = UnionSet(b1, bi)
        usa = UnionSetArray([b1, bi])

        # ------------------------------------------------------------------
        # plot using epsilon-close approximation (default threshold ε value)
        # ------------------------------------------------------------------
        plot(b2)
        plot(bp)
        plot(el)
        if n == 2
            plot(ncp)
        end
        @static if isdefined(@__MODULE__, :ExponentialUtilities) || isdefined(@__MODULE__, :Expokit)
            plot(em)
            plot(epm)
        end
        @static if isdefined(@__MODULE__, :Optim)
            if N == Float64 || n == 1  # Float32 in 2D requires promotion
                plot(its)
            else
                @test_broken plot(its) isa Vector{RecipesBase.RecipeData}
            end
        end
        if N == Float64 || n == 1  # Float32 in 2D requires promotion
            plot(itsa)
        else
            @test_broken plot(itsa) isa Vector{RecipesBase.RecipeData}
        end
        plot(hpg)
        plot(hpgo)
        plot(hpt)
        plot(hpt_empty)
        if n == 2
            plot(vpg)
        end
        plot(vpt)
        plot(hs)
        plot(hp)
        plot(hph)
        plot(uni)
        plot(l)
        if n == 2
            plot(l2D)
        end
        plot(sta)
        plot(sih)
        plot(lm)
        plot(rm)
        plot(ch)
        plot(cha)
        plot(ms)
        plot(msa)
        plot(cms)
        if n == 2
            plot(cp)
            plot(cpa)
        end

        @test_broken plot(dpz) isa Vector{RecipesBase.RecipeData}  # TODO DensePolynomialZonotope plot recipe is broken
        plot(spz)
        plot(sspz)
        plot(us)
        plot(usa)
    end

    # 3D
    @static if isdefined(@__MODULE__, :MiniQhull) && isdefined(@__MODULE__, :Polyhedra)
        n = 3
        bi = BallInf(zeros(N, 3), p1)
        t = Tetrahedron([N[0, 0, 0], N[1, 0, 0], N[0, 1, 0], N[0, 0, 1]])
        if N != Float32
            plot(bi)
        else
            # TODO the Float32 plot looks wrong (a pyramid) and randomly crashes in CI runs
            @test_broken false
        end
        plot(t)
    end
end

for N in [Float64]
    @static if isdefined(@__MODULE__, :StaticArrays)
        using StaticArrays: SA

        # test plot with static arrays input
        Z = Zonotope(SA[N(1), N(0)], SA[N(1) N(0); N(0) N(1)])
        plot(Z)
    end
end
