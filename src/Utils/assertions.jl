export activate_assertions,
       deactivate_assertions,
       are_assertions_enabled

# source: https://discourse.julialang.org/t/defensive-programming-assert/8383/11

# enable assertions by default
are_assertions_enabled() = true

# functions to (de)activate assertions
function activate_assertions()
    LazySets.eval(:(are_assertions_enabled() = true))
    nothing
end
function deactivate_assertions()
    LazySets.eval(:(are_assertions_enabled() = false))
    nothing
end

# override Base.@assert
macro assert(exs...)
    quote
        if $LazySets.are_assertions_enabled()
            Base.@assert $(map(esc, exs)...)
        end
    end
end
