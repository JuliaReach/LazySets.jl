using LazySets, Test
import Aqua, ExplicitImports

@testset "ExplicitImports tests" begin
    ignores = (:GLOBAL_RNG, :SamplerType, :invokelatest, :parse, :apply_recipe,
               :product, :checksquare, :Mesh, :coeff_table, :in_base, :numtype,
               :pos_table, :AbstractEnclosureAlgorithm, :_exp_remainder,
               :Symbolic, :Arrays, Symbol("@assert"))
    @test isnothing(ExplicitImports.check_all_explicit_imports_are_public(LazySets; ignore=ignores))
    @test isnothing(ExplicitImports.check_all_explicit_imports_via_owners(LazySets))
    ignores = (:Assertions, :Commutative, :Comparison, :EXACT, :Optimizer, :SIMPLEX, :Silent,
               :commutative, :uniontypes, :AbstractEnclosureAlgorithm, :EliminationAlgorithm,
               :Library, :get_degrees, :hcartesianproduct, :intersect, :isempty, :_exp_remainder,
               :hvectortype, :value, :setvrep!, :supportssolver, :vcartesianproduct, :Arr,
               :get_variables, :gradient, :Ellipsoid, :InteriorPoint, :Sets, :Translation,
               :ellipsoid, :invokelatest, :parse, :inf, :sup)
    @test isnothing(ExplicitImports.check_all_qualified_accesses_are_public(LazySets;
                                                                            ignore=ignores))
    ignores = (:inf, :sup)  # defined in IntervalArithmetic but imported through IntervalBoxes
    @test isnothing(ExplicitImports.check_all_qualified_accesses_via_owners(LazySets;
                                                                            ignore=ignores))
    @test isnothing(ExplicitImports.check_no_implicit_imports(LazySets))
    @test isnothing(ExplicitImports.check_no_self_qualified_accesses(LazySets))
    ignores = (:apply_recipe,)  # required for documentation
    @test isnothing(ExplicitImports.check_no_stale_explicit_imports(LazySets; ignore=ignores))
end

@testset "Aqua tests" begin
    Aqua.test_all(LazySets)
end
