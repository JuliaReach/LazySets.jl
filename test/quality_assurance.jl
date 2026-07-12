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
                                               :_σ_hyperplane_halfspace, :_σ_vertices, :_tohrep,
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
                                               :set_exponential_backend!,
                                               # RangeEnclosuresExt
                                               :_ρ_range_enclosures,
                                               # OptimExt
                                               :_line_search_optim,
                                               # MiniQhullExt
                                               :_triangulate_delaunay,
                                               # CDDLibExt
                                               :Library, :default_cddlib_backend,
                                               # IntervalBoxesExt
                                               :_difference,
                                               # GeometryBasicsExt
                                               :_area_polytope_3D, :_area_triangle_3D!,
                                               # IntervalMatricesExt
                                               :_exp_remainder, :taylor_expmap_remainder,
                                               # StaticArraysCoreExt
                                               :AbstractReductionMethod, :GIR05, :_convert_static,
                                               :_genmat_static, :_hcat_KLred, :_interval_hull,
                                               :_split_ret, :_to_colVector, :dir_east, :dir_north,
                                               :dir_south, :dir_west,
                                               # PolyhedraExt
                                               :EXACT, :EliminationAlgorithm, :GLPK_ON,
                                               :LinearMapElimination, :default_lp_solver_polyhedra,
                                               :default_polyhedra_backend_1d,
                                               :default_polyhedra_backend_nd, :hcartesianproduct,
                                               :hvectortype, :setvrep!, :supportssolver,
                                               :vcartesianproduct, :_area_Polyhedra,
                                               :_backend_solver_nd,
                                               :_cartesian_product_hrep_polyhedra,
                                               :_cartesian_product_vrep, :_convex_hull,
                                               :_get_elimination_instance,
                                               :_isempty_polyhedron_polyhedra,
                                               :_is_polyhedra_backend,
                                               :_minkowski_sum_hrep_preprocess,
                                               :_removehredundancy!, :_removevredundancy!,
                                               :_remove_redundant_vertices, :_vertices_list)
    ignores_all_explicit_imports_via_owners = (:BasicSymbolic,)
    ignores_all_qualified_accesses_are_public = (:Assertions, :Commutative, :Comparison, :EXACT,
                                                 :Optimizer, :SIMPLEX, :Silent, :commutative,
                                                 :uniontypes, :get_degrees, :intersect, :isempty,
                                                 :Ellipsoid, :inf, :sup, :OptimizerWithAttributes,
                                                 # fixed in versions after v"1.10"
                                                 :get_extension, :parse)
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
                                          no_stale_explicit_imports=(ignore=ignores_no_stale_explicit_imports,))
end

@safetestset "Aqua tests" begin
    import LazySets, Aqua

    Aqua.test_all(LazySets)
end
