const zero_itv_store = Dict{Type,IA.Interval}()
function zero_itv(N)
    res = get(zero_itv_store, N, nothing)
    if isnothing(res)
        res = IA.interval(zero(N))
        zero_itv_store[N] = res
    end
    return res
end

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

const sym_itv_store = Dict{Type,IA.Interval}()
function sym_itv(N)
    res = get(sym_itv_store, N, nothing)
    if isnothing(res)
        res = IA.interval(-one(N), one(N))
        sym_itv_store[N] = res
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

"""
    fast_interval_pow(a::IntervalArithmetic.Interval, n::Int)

Compute the `n`th power of an interval without using correct rounding.

### Input

- `a` -- interval (from `IntervalArithmetic.jl`)
- `n` -- integer

### Output

A non-rigorous approximation of `a^n`.

### Notes

For a rigorous approximation with correct rounding, use `a^n` from
`IntervalArithmetic.jl`.
"""
function fast_interval_pow(a::IA.Interval, n::Int)
    # TODO review after IntervalArithmetic.jl#388
    if iszero(n)
        return one(a)
    elseif isodd(n)
        return IA.interval(a.lo^n, a.hi^n)
    else
        if 0 ∈ a
            return IA.interval(zero(a.lo), max(abs(a.lo), abs(a.hi))^n)
        else
            lon = a.lo^n
            hin = a.hi^n
            return IA.interval(min(lon, hin), max(lon, hin))
        end
    end
end
