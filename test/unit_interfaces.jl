# check that all interace functions are provided

# --- LazySet ---

# support vector
@test check_method_implementation(LazySet, σ,
                                  Function[S -> (AbstractVector{Float64}, S)])
# dimension
@test check_method_implementation(LazySet, dim, Function[S -> (S,)])

# --- AbstractPolytope ---

# vertices list
@test check_method_implementation(AbstractPolytope, vertices_list,
                                  Function[S -> (S{Float64},)])
@test check_method_implementation(AbstractPointSymmetricPolytope, vertices_list,
                                  Function[S -> (S{Float64},)])

# --- AbstractPointSymmetric ---

# center
@test check_method_implementation(AbstractPointSymmetric, center,
                                  Function[S -> (S{Float64},)])
@test check_method_implementation(AbstractPointSymmetricPolytope, center,
                                  Function[S -> (S{Float64},)])

# --- AbstractHyperrectangle ---

# radius_hyperrectangle (2x)
@test check_method_implementation(AbstractHyperrectangle, radius_hyperrectangle,
                                  Function[S -> (S{Float64},)])
@test check_method_implementation(AbstractHyperrectangle, radius_hyperrectangle,
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
