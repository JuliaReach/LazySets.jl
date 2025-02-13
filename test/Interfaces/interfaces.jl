# check that all interface functions are implemented

# --- LazySet ---

# dimension
@test check_method_implementation(LazySet, dim, Function[S -> (S{Float64},)])

# concretize
@test check_method_implementation(LazySet, concretize, Function[S -> (S{Float64},)])

# --- ConvexSet ---

# support vector
@test check_method_implementation(ConvexSet, σ,
                                  Function[S -> (Vector{Float64}, S{Float64})])

# --- AbstractPolyhedron ---

# constraints list
@test check_method_implementation(AbstractPolyhedron, constraints_list,
                                  Function[S -> (S{Float64},)])

# --- AbstractPolytope ---

# vertices list
@test check_method_implementation(AbstractPolytope, vertices_list,
                                  Function[S -> (S{Float64},)])

# --- AbstractCentrallySymmetric and AbstractCentrallySymmetricPolytope ---

for T in (AbstractCentrallySymmetric, AbstractCentrallySymmetricPolytope)
    # center
    @test check_method_implementation(T, center,
                                      Function[S -> (S{Float64},)])
end

# --- AbstractZonotope ---

# genmat
@test check_method_implementation(AbstractZonotope, genmat,
                                  Function[S -> (S{Float64},)])
# generators
@test check_method_implementation(AbstractZonotope, generators,
                                  Function[S -> (S{Float64},)])

# --- AbstractHyperrectangle ---

# radius_hyperrectangle
@test check_method_implementation(AbstractHyperrectangle, radius_hyperrectangle,
                                  Function[S -> (S{Float64},)])

# --- AbstractSingleton ---

# element
@test check_method_implementation(AbstractSingleton, element,
                                  Function[S -> (S{Float64},)])

# --- AbstractPolygon ---

# tohrep
@test check_method_implementation(AbstractPolygon, tohrep,
                                  Function[S -> (S{Float64},)])
# tovrep
@test check_method_implementation(AbstractPolygon, tovrep,
                                  Function[S -> (S{Float64},)])

# --- AbstractHPolygon ---

# no functions yet

# --- AbstractAffineMap ---

# matrix
@test check_method_implementation(AbstractAffineMap, matrix,
                                  Function[S -> (S{Float64},)])
# vector
@test check_method_implementation(AbstractAffineMap, vector,
                                  Function[S -> (S{Float64},)])
# set
@test check_method_implementation(AbstractAffineMap, set,
                                  Function[S -> (S{Float64},)])

# --- AbstractBallp ---

# ball_norm
@test check_method_implementation(AbstractBallp, LazySets.ball_norm,
                                  Function[S -> (S{Float64},)])
# radius_ball
@test check_method_implementation(AbstractBallp, LazySets.radius_ball,
                                  Function[S -> (S{Float64},)])

# --- AbstractArraySet ---

for T in Base.uniontypes(LazySets.AbstractArraySet)
    # array
    @test check_method_implementation(T, array,
                                      Function[S -> (S{Float64},)])
end

# --- AbstractPolynomialZonotope ---

# center
@test check_method_implementation(AbstractPolynomialZonotope, center,
                                  Function[S -> (S{Float64},)])
# polynomial_order
@test check_method_implementation(AbstractPolynomialZonotope, polynomial_order,
                                  Function[S -> (S{Float64},)])
# ngens_dep
@test check_method_implementation(AbstractPolynomialZonotope, ngens_dep,
                                  Function[S -> (S{Float64},)])
# ngens_indep
@test check_method_implementation(AbstractPolynomialZonotope, ngens_indep,
                                  Function[S -> (S{Float64},)])

# --- AbstractSparsePolynomialZonotope ---

# expmat
@test check_method_implementation(AbstractSparsePolynomialZonotope, expmat,
                                  Function[S -> (S{Float64},)])
# genmat_dep
@test check_method_implementation(AbstractSparsePolynomialZonotope, genmat_dep,
                                  Function[S -> (S{Float64},)])
# genmat_indep
@test check_method_implementation(AbstractSparsePolynomialZonotope, genmat_indep,
                                  Function[S -> (S{Float64},)])
