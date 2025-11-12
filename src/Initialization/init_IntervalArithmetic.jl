const zero_itv_store = Dict{Type,IA.Interval}()
function zero_itv(N)
    res = get(zero_itv_store, N, nothing)
    if isnothing(res)
        res = IA.interval(zero(N))
        zero_itv_store[N] = res
    end
    return res
end

const sym_itv_store = Dict{Type,IA.Interval}()
function sym_itv(N)
    res = get(sym_itv_store, N, nothing)
    if isnothing(res)
        res = IA.interval(-one(N), one(N))
        sym_itv_store[N] = res
    end
    return res
end
