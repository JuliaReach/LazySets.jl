# check that the following functions are provided by all LazySet types
for S in subtypes(LazySet)
    # support vector
    @test method_exists(σ, (AbstractVector{Float64}, S)) ||
          method_exists(σ, (AbstractVector{Float64}, S{Float64}))
    # dimension
    @test method_exists(dim, (S,))
end
