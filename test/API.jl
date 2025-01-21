@testset "API loading" begin
    using LazySets.API
end

@testset "API set" begin
    struct MySet <: API.LazySet end
    MySet()
end

@testset "API operations" begin
    X = MySet()

    # unary set operations
    for f in (an_element, area, center, complement, concretize, constraints_list, constraints,
              convex_hull, diameter, dim, eltype, extrema, high, ispolyhedral, isbounded, isempty,
              isoperation, isuniversal, low, norm, radius, rectify, reflect, sample, surface,
              vertices_list, vertices, volume)
        @test isnothing(f(X))
    end

    # unary set operations with optional arguments
    for f in (isempty, isuniversal)
        @test isnothing(f(X, true))
    end
    for f in (diameter, norm, radius)
        @test isnothing(f(X, 1))
    end

    # unary set-type operations
    for f in (eltype, isboundedtype, isconvextype, isoperationtype, rand)
        @test isnothing(f(MySet))
    end

    # mixed set operations (typically with vectors or matrices)
    n = 1
    for f in (scale!, scale)
        @test isnothing(f(n, X))
    end
    for f in (center, extrema, high, low, sample)
        @test isnothing(f(X, n))
    end
    v = [1]
    for f in (distance, permute, project, translate!, translate)
        @test isnothing(f(X, v))
    end
    for f in (distance, is_interior_point, ∈, support_function, ρ, support_vector, σ)
        @test isnothing(f(v, X))
    end
    M = hcat(1)
    for f in (exponential_map, linear_map)
        @test isnothing(f(M, X))
    end
    @test isnothing(affine_map(M, X, v))

    # binary set operations
    for f in (cartesian_product, difference, distance, exact_sum, intersection,
              is_intersection_empty, isdisjoint, ≈, ==, ⊆, isequivalent, ⊂, linear_combination,
              minkowski_difference, pontryagin_difference, minkowski_sum)
        @test isnothing(f(X, X))
    end
    @test isnothing(distance(X, X; p=2.0))
end
