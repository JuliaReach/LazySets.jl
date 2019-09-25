export activate_benchmark_mode,
       deactivate_benchmark_mode

# source: https://discourse.julialang.org/t/defensive-programming-assert/8383/11

# functions to (de)activate assertions
function activate_benchmark_mode()
    LazySets.eval(:(_debug_disabled() = true))
    nothing
end
function deactivate_benchmark_mode()
    LazySets.eval(:(_debug_disabled() = false))
    nothing
end

# enable assertions by default
deactivate_benchmark_mode()

# override Base.@assert
macro assert(exs...)
    quote
        if !$LazySets._debug_disabled()
            Base.@assert $(map(esc, exs)...)
        end
    end
end
