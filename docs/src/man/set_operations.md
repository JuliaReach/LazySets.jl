# Operations on sets

In this section we show which typical set operations this library supports.

```@contents
Pages = ["set_operations.md"]
Depth = 3
```

We use the following four sets for illustration.

```@example set_operations
using LazySets, LazySets.Approximations, Plots
B1 = Ball1(-ones(2), 1.)
B2 = Ball2(ones(2), 1.)
BI = BallInf(zeros(2), 1.)
H = Hyperrectangle(ones(2), ones(2))
sets = [B1, B2, BI, H]

function plot_sets(sets)
    for S in sets
        println(S)
        plot!(S, 1e-2, fillalpha=0.1)
    end
end

function plot_points(points, prefix)
    for i in eachindex(points)
        p = points[i]
        num_occur = length(find(x -> x == p, points[1:i]))
        x = p[1]
        y = p[2]
        if num_occur == 1
            x += 0.15
        elseif num_occur == 2
            y += 0.15
        elseif num_occur == 3
            x -= 0.15
        else
            y -= 0.15
        end
        plot!(Singleton(p))
        plot!(annotations=(x, y, text("$(prefix)$(i)")))
    end
end

plot1 = plot()
plot_sets(sets)
plot1
```

## Unary operations

The following table lists all operations (apart from `dim` and `σ`, which must
be supported by all types) that take one convex set as argument in the columns.
In the rows we list all set types, both the interfaces (where we abbreviate the
`Abstract` prefix), the basic set types, and the lazy set operations, each
sorted alphabetically.
The table entries have the following meaning.
- An "x" indicates that the operation is implemented for the respective set.
- "(⋅)" indicates that the operation is inherited from a supertype.


| type ↓ \ operation →         | radius | diameter | norm | an_element | ∈ |
|------------------------------|--------|----------|------|------------|---|
| **Interfaces**               |        |          |      |            |   |
| `LazySet`                    |        |          |      |            |   |
| `AHPolygon`                  |        |          |      |     x      | x |
| `AHyperrectangle`            |   x    |    x     |  x   |    (x)     | x |
| `APointSymmetric`            |        |          |      |     x      |   |
| `APointSymmetricPolytope`    |        |          |      |     x      |   |
| `APolygon`                   |        |          |      |            |   |
| `APolytope`                  |        |          |      |            |   |
| `ASingleton`                 |  (x)   |   (x)    | (x)  |     x      | x |
|                              |        |          |      |            |   |
| **Basic set types**          |        |          |      |            |   |
| `Ball1`                      |        |          |      |    (x)     | x |
| `Ball2`                      |        |          |      |    (x)     | x |
| `BallInf`                    |   x    |   (x)    | (x)  |    (x)     |(x)|
| `Ballp`                      |        |          |      |    (x)     | x |
| `Ellipsoid`                  |        |          |      |    (x)     | x |
| `EmptySet`                   |        |          |      |     x      | x |
| `HPolygon`/`HPolygonOpt`     |        |          |      |    (x)     |(x)|
| `HPolytope`                  |        |          |      |            | x |
| `Hyperrectangle`             |  (x)   |   (x)    | (x)  |    (x)     |(x)|
| `Singleton`                  |  (x)   |   (x)    | (x)  |    (x)     |(x)|
| `VPolygon`                   |        |          |      |     x      | x |
| `ZeroSet`                    |  (x)   |   (x)    | (x)  |    (x)     | x |
| `Zonotope`                   |        |          |      |    (x)     | x |
|                              |        |          |      |            |   |
| **Lazy set operation types** |        |          |      |            |   |
| `CartesianProduct`           |        |          |      |            | x |
| `CartesianProductArray`      |        |          |      |            | x |
| `ConvexHull`                 |        |          |      |            |   |
| `ConvexHullArray`            |        |          |      |            |   |
| `ExponentialMap`             |        |          |      |            | x |
| `ExponentialProjectionMap`   |        |          |      |            |   |
| `HalfSpace`                  |        |          |      |     x      | x |
| `Hyperplane`                 |        |          |      |     x      | x |
| `Intersection`               |        |          |      |            | x |
| `LinearMap`                  |        |          |      |            | x |
| `MinkowskiSum`               |        |          |      |            |   |
| `MinkowskiSumArray`          |        |          |      |            |   |
| `SymmetricIntervalHull`      |  (x)   |   (x)    | (x)  |    (x)     |(x)|


### `radius`

This function returns the radius of a set.
It is defined as the radius of the enclosing ball (of the given norm) of
minimal volume with the same center.

```@example set_operations
radius(B1), radius(B2), radius(BI), radius(H)
```

### `diameter`

This function returns the diameter of a set.
It is defined as the diameter of the enclosing ball (of the given norm) of
minimal volume with the same center.
The implementation is inherited for all set types if the norm is the infinity
norm, in which case the result is defined as twice the radius.

```@example set_operations
diameter(B1), diameter(B2), diameter(BI), diameter(H)
```

### `norm`

This function returns the norm of a set.
It is defined as the norm of the enclosing ball (of the given norm) of minimal
volume centered in the origin.

```@example set_operations
# print 1-norm, 2-norm, and infinity norm (if available)
println(("-", "-", norm(B1, Inf)))
println(("-", "-", norm(B2, Inf)))
println((norm(BI, 1), norm(BI, 2), norm(BI, Inf)))
println((norm(H, 1), norm(H, 2), norm(H, Inf)))
```

### `an_element`

This function returns some element in the set.
Consecutive calls to this function typically return the same element.

```@example set_operations
an_element(B1), an_element(B2), an_element(BI), an_element(H)
```

### `∈`

This function checks containment of a given vector in the set.
The operator can be used in infix notation (`v ∈ S`) and in inverse operand
order (`S ∋ v`).
Alias: `in`

```@example set_operations
p1 = [1.5, 1.5]
p2 = [0.1, 0.1]
p3 = [-0.9, -0.8]
points = [p1, p2, p3]

for p in [p1, p2, p3]
    println("$p ∈ (B1, B2, BI, H)? ($(p ∈ B1), $(p ∈ B2), $(p ∈ BI), $(p ∈ H))")
end
```

```@example set_operations
plot1 = plot()
plot_sets(sets)
plot_points(points, "p")
plot1
```


## Binary operations

The following table lists all operations that take two convex set as argument in
the entries.
In the rows we list all set types, both the interfaces (where we abbreviate the
`Abstract` prefix), the basic set types, and the lazy set operations, each
sorted alphabetically.
In the columns we also list the operations, but abbreviated.
The table entries consist of subsets of the following list of operations.
- "⊆" stands for the subset check operation.
- "∩=∅" stands for the intersection emptiness check operation.
- "(⋅)" indicates that the operation is inherited from a supertype.
- "[⋅]" indicates that a possible type ambiguity can occur but is resolved.
- "{⋅}" indicates that a more efficient implementation is used instead of the
  inherited one.


| type ↓ \ type →              |LazyS|AHPon|AHrec|APSym|APSPol|APgon|APtop|ASingle| B1  | B2  |BInf | Bp  |Ellip|EmpSt|HPgon|HPtop|Hrect |Single |VPgon|ZeroSet|Zonot| CP  | CPA | CH  | CHA |EMap | EPM |HalfS|HypPl|Inter|LMap |MSum | MSA | SIH |
|------------------------------|-----|-----|-----|-----|------|-----|-----|-------|-----|-----|-----|-----|-----|-----|-----|-----|------|-------|-----|-------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| **Interfaces**               |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `LazySet`                    |     |     |  ⊆  |     |      |     |     |(⊆),∩=∅|     |     | (⊆) |     |     |     |     |     | (⊆)  | (⊆)   |     | (⊆)   |     |     |     |     |     |     |     |     |     |     |     |     |     | (⊆) |
| `AHPolygon`                  | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `AHyperrectangle`            | (⊆) | (⊆) |{⊆},∩=∅|(⊆)| (⊆)  | (⊆) | (⊆) |[⊆,∩=∅]| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |(⊆,∩=∅)|(⊆,∩=∅)|(⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `APointSymmetric`            |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `APointSymmetricPolytope`    | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `APolygon`                   | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `APolytope`                  |  ⊆  | (⊆) | [⊆] | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `ASingleton`                 |{⊆},∩=∅|(⊆,∩=∅)|[⊆,∩=∅]|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|[⊆,∩=∅]|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|((⊆,∩=∅))|(⊆,∩=∅)|((⊆,∩=∅))|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|
|                              |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| **Basic set types**          |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `Ball1`                      | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `Ball2`                      |     |     | (⊆) |     |      |     |     |{⊆},(∩=∅)|   |⊆,∩=∅|     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `BallInf`                    | (⊆) | (⊆) |(⊆,∩=∅)|(⊆)| (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `Ballp`                      |     |     | (⊆) |     |      |     |     |{⊆},(∩=∅)|   |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `Ellipsoid`                  |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `EmptySet`                   |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `HPolygon`/`HPolygonOpt`     | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `HPolytope`                  | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `Hyperrectangle`             | (⊆) | (⊆) |(⊆,∩=∅)|(⊆)| (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `Singleton`                  |(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|
| `VPolygon`                   | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `ZeroSet`                    |(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|
| `Zonotope`                   | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
|                              |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| **Lazy set operation types** |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `CartesianProduct`           |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `CartesianProductArray`      |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `ConvexHull`                 |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `ConvexHullArray`            |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `ExponentialMap`             |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `ExponentialProjectionMap`   |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `HalfSpace`                  |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `Hyperplane`                 |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `Intersection`               |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `LinearMap`                  |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `MinkowskiSum`               |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `MinkowskiSumArray`          |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `SymmetricIntervalHull`      | (⊆) | (⊆) |(⊆,∩=∅)|(⊆)| (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |


### `⊆`

This function checks whether a set is a subset of another set.
It can optionally produce a witness if the subset relation does not hold.
The operator can be used in infix notation (`X ⊆ S`).
Alias: `issubset`

```@example set_operations
println(B1 ⊆ B2)
w1 = ⊆(B1, B2, true)[2]
println(B1 ⊆ BI)
w2 = ⊆(B1, BI, true)[2]
println(B1 ⊆ H)
w3 = ⊆(B1, H, true)[2]
# 'B2 ⊆ B1' is not supported yet
# w11 = ⊆(B2, B1, true)[2]
println(B2 ⊆ BI)
w4 = ⊆(B2, BI, true)[2]
println(B2 ⊆ H)
println(BI ⊆ B1)
w5 = ⊆(BI, B1, true)[2]
println(BI ⊆ B2)
w6 = ⊆(BI, B2, true)[2]
println(BI ⊆ H)
w7 = ⊆(BI, H, true)[2]
println(H ⊆ B1)
w8 = ⊆(H, B1, true)[2]
println(H ⊆ B2)
w9 = ⊆(H, B2, true)[2]
println(H ⊆ BI)
w10 = ⊆(H, BI, true)[2];
```

```@example set_operations
witnesses = [w1, w2, w3, w4, w5, w6, w7, w8, w9, w10]

plot1 = plot()
plot_sets(sets)
plot_points(witnesses, "w")
plot1
```

### `is_intersection_empty`

This function checks whether the intersection of two sets is empty.
It can optionally produce a witness if the intersection is nonempty.

```@example set_operations
println(is_intersection_empty(BI, H))
w1 = is_intersection_empty(BI, H, true)[2]
# none of the other combinations are supported yet
# is_intersection_empty(B1, B2)
# is_intersection_empty(B1, BI)
# is_intersection_empty(B1, H)
# w2 = is_intersection_empty(B1, H, true)[2]
# is_intersection_empty(B2, BI)
# is_intersection_empty(B2, H)
```

```@example set_operations
witnesses = [w1]

plot1 = plot()
plot_sets(sets)
plot_points(witnesses, "w")
plot1
```
