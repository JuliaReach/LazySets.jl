using SafeTestsets

@safetestset "ExplicitImports tests" begin
    using Test
    import LazySets, ExplicitImports

    ignores_all_explicit_imports_are_public = (:GLOBAL_RNG, :SamplerType, :invokelatest, :parse,
                                               :apply_recipe, :product, :checksquare, :Mesh,
                                               :coeff_table, :in_base, :numtype, :pos_table,
                                               :AbstractEnclosureAlgorithm, :_exp_remainder,
                                               :BasicSymbolic, :Symbolic, :Arrays,
                                               Symbol("@assert"))
    ignores_all_explicit_imports_via_owners = (:BasicSymbolic,)
    ignores_all_qualified_accesses_are_public = (:Assertions, :Commutative, :Comparison, :EXACT,
                                                 :Optimizer, :SIMPLEX, :Silent, :commutative,
                                                 :uniontypes, :AbstractEnclosureAlgorithm,
                                                 :EliminationAlgorithm, :Library, :get_degrees,
                                                 :hcartesianproduct, :intersect, :isempty,
                                                 :_exp_remainder, :hvectortype, :value, :setvrep!,
                                                 :supportssolver, :vcartesianproduct, :Arr,
                                                 :get_variables, :gradient, :Ellipsoid,
                                                 :InteriorPoint, :Sets, :Translation, :ellipsoid,
                                                 :invokelatest, :parse, :inf, :sup)
    ignores_all_qualified_accesses_via_owners = (:inf, :sup)  # defined in IntervalArithmetic but imported through IntervalBoxes
    ignores_no_stale_explicit_imports = (:apply_recipe,)  # required for documentation
    ExplicitImports.test_explicit_imports(LazySets;
                                          all_explicit_imports_are_public=(ignore=ignores_all_explicit_imports_are_public,),
                                          all_explicit_imports_via_owners=(ignore=ignores_all_explicit_imports_via_owners,),
                                          all_qualified_accesses_are_public=(ignore=ignores_all_qualified_accesses_are_public,),
                                          all_qualified_accesses_via_owners=(ignore=ignores_all_qualified_accesses_via_owners,),
                                          no_stale_explicit_imports=(ignore=ignores_no_stale_explicit_imports,))
end

@safetestset "Aqua tests" begin
    import LazySets, Aqua

    Aqua.test_all(LazySets)
end
