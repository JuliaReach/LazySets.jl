"""
    read_gen(filename::String)

Read a sequence of polygons stored in vertex representation (gen format).

### Input

- `filename` -- path of the file containing the polygons

### Output

A list of polygons in vertex representation.

### Notes

The `x` and `y` coordinates of each vertex should be separated by an empty space
and polygons are separated by empty lines (even the last polygon). For example:

```julia
1.01 1.01
0.99 1.01
0.99 0.99
1.01 0.99

0.908463 1.31047
0.873089 1.31047
0.873089 1.28452
0.908463 1.28452


```
This is parsed as

```julia
2-element Array{VPolygon{Float64, Vector{Float64}},1}:
 VPolygon{Float64, Vector{Float64}}([[1.01, 1.01], [0.99, 1.01], [0.99, 0.99], [1.01, 0.99]])
 VPolygon{Float64, Vector{Float64}}([[0.908463, 1.31047], [0.873089, 1.31047], [0.873089, 1.28452], [0.908463, 1.28452]])
```
"""
function read_gen(filename::String)
    Mi = Vector{Vector{Float64}}()
    P = Vector{VPolygon{Float64, Vector{Float64}}}()

    # detect when we finished reading a new polygon; needed because polygons may
    # be separated by more than one end-of-line
    new_polygon = true
    open(filename) do f
        for line in eachline(f)
            if !isempty(line)
                push!(Mi, map(x -> parse(Float64, x), split(line)))
                new_polygon = true
             elseif isempty(line) && new_polygon
                push!(P, VPolygon(Mi))
                Mi = Vector{Vector{Float64}}()
                new_polygon = false
             end
        end
    end
    return P
end
