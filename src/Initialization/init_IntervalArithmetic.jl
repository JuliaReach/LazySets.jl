const zero_itv_store = Dict{Type,IA.Interval}()
function zero_itv(N)
    res = get(zero_itv_store, N, nothing)
    if isnothing(res)
        res = IA.interval(zero(N))
        zero_itv_store[N] = res
    end
    return res
end

const zero_box_store = Dict{Type,Dict{Integer,Vector{<:IA.Interval}}}()
function zero_box(n::Int, N=Float64)
    d = get(zero_box_store, N, nothing)
    if isnothing(d)
        d = Dict{Integer,Vector{<:IA.Interval}}()
        zero_box_store[N] = d
    end
    res = get(d, n, nothing)
    if isnothing(res)
        res = fill(zero_itv(N), n)
        d[n] = res
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

const sym_box_store = Dict{Type,Dict{Integer,Vector{<:IA.Interval}}}()
function sym_box(n::Int, N=Float64)
    d = get(sym_box_store, N, nothing)
    if isnothing(d)
        d = Dict{Integer,Vector{<:IA.Interval}}()
        sym_box_store[N] = d
    end
    res = get(d, n, nothing)
    if isnothing(res)
        res = fill(sym_itv(N), n)
        d[n] = res
    end
    return res
end
