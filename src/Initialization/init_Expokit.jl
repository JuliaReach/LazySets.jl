eval(quote
    using .Expokit: expmv
end)

eval(load_expokit_sparsematrixexp())
eval(load_expokit_exponentialmap())
eval(load_expokit_exponentialprojectionmap())
