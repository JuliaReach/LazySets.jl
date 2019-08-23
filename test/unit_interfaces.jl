# check that all interace functions are provided

# --- LazySet ---

# support vector
@test check_method_implementation(LazySet, σ,
                                  Function[S -> (Vector{Float64}, S{Float64})])
# dimension
check_method_implementation(LazySet, dim, Function[S -> (S{Float64},)])

# --- AbstractPolyhedron ---

# constraints list
@test check_method_implementation(AbstractPolyhedron, constraints_list,
                                  Function[S -> (S{Float64},)])

# --- AbstractPolytope ---

# vertices list
@test check_method_implementation(AbstractPolytope, vertices_list,
                                  Function[S -> (S{Float64},)])

# --- AbstractCentrallySymmetric ---

# center
@test check_method_implementation(AbstractCentrallySymmetric, center,
                                  Function[S -> (S{Float64},)])

# --- AbstractCentrallySymmetricPolytope ---
# (copies methods from AbstractCentrallySymmetric)

# center (from AbstractCentrallySymmetric)
@test check_method_implementation(AbstractCentrallySymmetricPolytope, center,
                                  Function[S -> (S{Float64},)])

# --- AbstractZonotope ---

# genmat
@test check_method_implementation(AbstractZonotope, genmat,
                                  Function[S -> (S{Float64},)])
# generators
@test check_method_implementation(AbstractZonotope, generators,
                                  Function[S -> (S{Float64},)])

# --- AbstractHyperrectangle ---

# radius_hyperrectangle (2x)
@test check_method_implementation(AbstractHyperrectangle, radius_hyperrectangle,
                                  Function[S -> (S{Float64},)])
@test check_method_implementation(AbstractHyperrectangle, radius_hyperrectangle,
                                  Function[S -> (S{Float64}, Int)])
# isflat
@test check_method_implementation(AbstractHyperrectangle, isflat,
                                  Function[S -> (S{Float64},)])

# --- AbstractSingleton ---

# element (2x)
@test check_method_implementation(AbstractSingleton, element,
                                  Function[S -> (S{Float64},)])
@test check_method_implementation(AbstractSingleton, element,
                                  Function[S -> (S{Float64}, Int)])

# --- AbstractPolygon ---

# tohrep
@test check_method_implementation(AbstractPolygon, tohrep,
                                  Function[S -> (S{Float64},)])
# tovrep
@test check_method_implementation(AbstractPolygon, tovrep,
                                  Function[S -> (S{Float64},)])

# --- AbstractHPolygon ---

# no functions yet
