const zero_box_store = Dict{Type,Dict{Integer,IA.IntervalBox}}()
function zero_box(n::Int, N=Float64)
    d = get(zero_box_store, N, nothing)
    if isnothing(d)
        d = Dict{Integer,IA.IntervalBox}()
        zero_box_store[N] = d
    end
    res = get(d, n, nothing)
    if isnothing(res)
        res = IA.IntervalBox(zero_itv(N), n)
        d[n] = res
    end
    return res
end

const sym_box_store = Dict{Type,Dict{Integer,IA.IntervalBox}}()
function sym_box(n::Int, N=Float64)
    d = get(sym_box_store, N, nothing)
    if isnothing(d)
        d = Dict{Integer,IA.IntervalBox}()
        sym_box_store[N] = d
    end
    res = get(d, n, nothing)
    if isnothing(res)
        res = IA.IntervalBox(sym_itv(N), n)
        d[n] = res
    end
    return res
end
