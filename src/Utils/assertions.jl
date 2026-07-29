# convenience functions to (de)activate assertions
function activate_assertions()
    for m in [LazySets, Approximations]
        ReachabilityBase.Assertions.activate_assertions(m)
    end
end

function deactivate_assertions()
    for m in [LazySets, Approximations]
        ReachabilityBase.Assertions.deactivate_assertions(m)
    end
end
