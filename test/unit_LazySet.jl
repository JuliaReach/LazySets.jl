# check that the following functions are provided by all LazySet types
for S in subtypes(LazySet)
    # support vector
    @test method_exists(σ, (Vector{Float64}, S))
    # dimension
    @test method_exists(dim, (S,))
end
