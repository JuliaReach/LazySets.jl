using LazySets, Test
import Aqua

@testset "Aqua tests" begin
    Aqua.test_all(LazySets; ambiguities=false)

    # do not warn about ambiguities in dependencies
    Aqua.test_ambiguities(LazySets)
end
