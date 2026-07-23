using SafeTestsets

@safetestset "ExplicitImports tests" begin
    using Test
    import LazySets, ExplicitImports

    ignores_all_explicit_imports_are_public = (:GLOBAL_RNG, :SamplerType, :invokelatest, :parse,
                                               :apply_recipe, :product, :checksquare, :Mesh,
                                               :coeff_table, :in_base, :numtype, :pos_table,
                                               :AbstractEnclosureAlgorithm, :_exp_remainder,
                                               :BasicSymbolic, :Symbolic, :Arrays,
                                               Symbol("@assert"),
                                               # package extensions
                                               # LazySetsExt
                                               :AbstractLinearMapAlgorithm, :HPoly, :LinearMapVRep,
                                               :STAR, :default_polyhedra_backend,
                                               :_Ballp_special_cases, :_HPolyhedron,
                                               :_absorb_generators_spz, :_an_element_halfspace,
                                               :_an_element_helper_hyperplane, :_center,
                                               :_constraints_list_Vector,
                                               :_constraints_list_hparallelotope,
                                               :_constraints_list_hyperplane,
                                               :_constraints_list_singleton_Vector,
                                               :_convert_HPolygon, :_difference_universe2,
                                               :_infeasible_constraints_list,
                                               :_intersection_vrep_2d,
                                               :_isdisjoint_hyperplane_hyperplane,
                                               :_linear_map_hrep, :_linear_map_hrep_helper,
                                               :_linear_map_vrep, :_linear_map_zonotope,
                                               :_minkowski_difference_universe2,
                                               :_non_element_halfspace, :_normalize_halfspace,
                                               :_sort_constraints, :_ρ_vertices,
                                               :_σ_hyperplane_halfspace, :_σ_vertices,
                                               :_witness_result_empty, Symbol("@validate"),
                                               Symbol("@validate_commutative"),
                                               # RecipesBaseExt
                                               :plot_recipe, :plot_vlist, :_plot_recipe_3d_polytope,
                                               # DistributionsExt
                                               :RejectionSampler, :_sample_unit_nball_muller!,
                                               :_sample_unit_nsphere_muller!,
                                               # SCSExt
                                               :MOI, :Optimizer, :sdp_solver, :set_sdp_solver!,
                                               :_default_sdp_solver,
                                               # IpoptExt
                                               :_default_nln_solver,
                                               # MakieExt
                                               :Automatic,
                                               # IntervalConstraintProgrammingExt
                                               :_contract_zonotope_halfspace_ICP,
                                               # SetProgExt
                                               :InteriorPoint, :Translation, :ellipsoid,
                                               :_underapproximate_ellipsoid, :default_sdp_solver,
                                               # SymbolicsExt (not needed in v"1.12")
                                               :Arr, :get_variables, :gradient, :value,
                                               # ExpokitExt
                                               :_expmv, :exponential_backend,
                                               :set_exponential_backend!)
    ignores_all_explicit_imports_via_owners = (:BasicSymbolic,)
    ignores_all_qualified_accesses_are_public = (:Assertions, :Commutative, :Comparison, :EXACT,
                                                 :Optimizer, :SIMPLEX, :Silent, :commutative,
                                                 :uniontypes, :AbstractEnclosureAlgorithm,
                                                 :EliminationAlgorithm, :Library, :get_degrees,
                                                 :hcartesianproduct, :intersect, :isempty,
                                                 :_exp_remainder, :hvectortype, :value, :setvrep!,
                                                 :supportssolver, :vcartesianproduct,
                                                 :get_variables, :gradient, :Ellipsoid,
                                                 :invokelatest, :parse, :inf, :sup,
                                                 :OptimizerWithAttributes,
                                                 # fixed in versions after v"1.10"
                                                 :get_extension)
    ignores_all_qualified_accesses_via_owners = (:inf, :sup)  # defined in IntervalArithmetic but imported through IntervalBoxes
    ignores_no_stale_explicit_imports = (
                                         # required for documentation
                                         :apply_recipe,
                                         # IntervalConstraintProgrammingExt
                                         Symbol("@variables")
                                         )
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
