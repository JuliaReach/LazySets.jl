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
        num_occur = length(findfirst(x -> x == p, points[1:i]))
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

The following table lists all operations that take one convex set as argument in
the columns.
In the rows we list all set types, both the interfaces (where we abbreviate the
`Abstract` prefix), the basic set types, and the lazy set operations, each
sorted alphabetically.
The table entries have the following meaning.
- "x" indicates that the operation is implemented for the respective set type.
- "i" indicates that the operation is inherited from a supertype.
- "(·)" indicates that the operation is partly implemented/inherited.


| type ↓ \ operation →         | dim | ρ | σ | an_element | ∈ | isempty | isbounded | linear_map | translate | norm | radius | diameter |
|------------------------------|-----|---|---|------------|---|---------|-----------|------------|-----------|------|--------|----------|
| **Interfaces**               |     |   |   |            |   |         |           |            |           |      |        |          |
| `LazySet`                    |     | x |   | x          |   |         | x         |            |           |      |        | x        |
| `APolytope`                  |     | i |   | i          |   | x       | x         | x          |           |      |        | i        |
| `ACentrallySymmetric`        | x   | i |   | x          |   | x       | x         |            |           |      |        | i        |
| `ACentrallySymmetricPolytope`| i   | i |   | i          |   | x       | i         | i          |           |      |        | i        |
| `APolygon`                   | x   | i |   | i          |   | i       | i         | i          |           |      |        | i        |
| `AHyperrectangle`            | i   | i | x | i          | x | i       | i         | i          |           | x    | x      | i        |
| `AHPolygon`                  | i   | i |   | x          | x | i       | i         | i          |           |      |        | i        |
| `ASingleton`                 | i   | i | x | i          | x | i       | i         | x          |           | i    | i      | i        |
|                              |     |   |   |            |   |         |           |            |           |      |        |          |
| **Basic set types**          |     |   |   |            |   |         |           |            |           |      |        |          |
| `Ball1`                      | i   | i | x | i          | x | i       | i         | i          | x         |      |        | i        |
| `Ball2`                      | i   | i | x | i          | x | i       | i         |            | x         |      |        | i        |
| `BallInf`                    | i   | i | i | i          | i | i       | i         | i          | x         | i    | x      | i        |
| `Ballp`                      | i   | i | x | i          | x | i       | i         |            | x         |      |        | i        |
| `Ellipsoid`                  | i   | x | x | i          | x | i       | i         |            | x         |      |        | i        |
| `EmptySet`                   | x   | i | x | x          | x | x       | x         |            | x         | x    | x      | x        |
| `HalfSpace`                  | x   | x | x | x          | x | x       | x         |            | x         |      |        | i        |
| `HPolygon`/`HPolygonOpt`     | i   | i | x | i          | i | i       | i         | i          | x         |      |        | i        |
| `HPolyhedron`                | x   | x | x | i          | x | x       | x         | x          | x         |      |        | i        |
| `HPolytope`                  | x   | x | x | i          | x | x       | i         | x          | x         |      |        | i        |
| `Hyperplane`                 | x   | x | x | x          | x | x       | x         |            | x         |      |        | i        |
| `Hyperrectangle`             | i   | i | i | i          | i | i       | i         | i          | x         | i    | i      | i        |
| `Interval`                   | x   | i | x | x          | x | i       | i         | i          | x         | i    | i      | i        |
| `Line`                       | x   | i | x | x          | x | x       | x         |            | x         |      |        | i        |
| `LineSegment`                | x   | i | x | x          | x | i       | i         | i          | x         |      |        | i        |
| `Singleton`                  | i   | i | i | i          | i | i       | i         | i          | x         | i    | i      | i        |
| `Universe`                   | x   | x | x | x          | x | x       | x         |            | x         | x    | x      | x        |
| `VPolygon`                   | i   | i | x | x          | x | i       | i         | x          | x         |      |        | i        |
| `VPolytope`                  | x   | i | x | i          | x | i       | i         | x          | x         |      |        | i        |
| `ZeroSet`                    | x   | i | x | i          | x | i       | i         | x          | x         | i    | i      | i        |
| `Zonotope`                   | i   | x | x | i          | x | i       | i         | x          | x         |      |        | i        |
|                              |     |   |   |            |   |         |           |            |           |      |        |          |
| **Lazy set operation types** |     |   |   |            |   |         |           |            |           |      |        |          |
| `CartesianProduct`           | x   | x | x | i          | x | x       | x         |            |           |      |        | i        |
| `CartesianProductArray`      | x   | x | x | i          | x | x       | x         |            |           |      |        | i        |
| `ConvexHull`                 | x   | x | x | i          |   | x       | x         |            |           |      |        | i        |
| `ConvexHullArray`            | x   | x | x | i          |   | x       | x         |            |           |      |        | i        |
| `ExponentialMap`             | x   | x | x | i          | x | x       | x         |            |           |      |        | i        |
| `ExponentialProjectionMap`   | x   | i | x | i          |   | x       | x         |            |           |      |        | i        |
| `Intersection`               | x   | x |   | i          | x | x       | x         |            |           |      |        | i        |
| `IntersectionArray`          | x   | i |   | i          | x |         | x         |            |           |      |        | i        |
| `LinearMap`                  | x   | x | x | x          | x | x       | x         |            |           |      |        | i        |
| `MinkowskiSum`               | x   | x | x | i          |   | x       | x         |            |           |      |        | i        |
| `MinkowskiSumArray`          | x   | x | x | i          |   | x       | x         |            |           |      |        | i        |
| `CacheMinkowskiSum`          | x   | i | x | i          |   | x       | x         |            |           |      |        | i        |
| `ResetMap`                   | x   | x | x | x          |   | x       |           |            |           |      |        | i        |
| `SymmetricIntervalHull`      | x   | i | x | i          | i | i       | i         | i          |           | i    | i      | i        |
|                              |     |   |   |            |   |         |           |            |           |      |        |          |
| **Non-convex operations**    |     |   |   |            |   |         |           |            |           |      |        |          |
| `Complement`                 | x   |   |   |            | x | x       |           |            |           |      |        |          |
| `Rectification`              | x   | i |(x)|            |   |         |           |            |           |      |        |          |
| `UnionSet`                   | x   | x | x | x          | x | x       | x         |            |           |      |        |          |
| `UnionSetArray`              | x   | x | x | x          | x | x       | x         |            |           |      |        |          |


### `dim`

This function returns the dimension of the set.

```@example set_operations
dim(B1), dim(B2), dim(BI), dim(H)
```

### `ρ`/`σ`

These functions return the support function resp. the support vector of the set.

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

### `isempty`

This function checks if the set is empty.

### `linear_map`

This function applies a concrete linear map to the set.
The resulting set may be of a different type.

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


## Binary operations

The following table lists all operations that take two convex set as argument in
the entries.
In the rows we list all set types, both the interfaces (where we abbreviate the
`Abstract` prefix), the basic set types, and the lazy set operations, each
sorted alphabetically.
In the columns we also list the operations, but abbreviated.
The table entries consist of subsets of the following list of operations.
- "⊆" stands for the subset check `issubset`.
- "⊎" stands for the disjointness check `isdisjoint`.
- "∩" stands for the concrete intersection operation `intersection`.
- "C" stands for the conversion operation `convert`.
- "-" indicates that the two types' dimensionality constraints are incompatible.
- A suffix "i" indicates that the operation is inherited from a supertype.


| type ↓ \ type →               |LazyS      |APtop      |ACSym      |ACSPt      |APgon      |AHrec      |AHPgn      |ASing      |Ball1      |Ball2      |BInf       |Ballp      |Ellip      |Empty      |HalfS      |HPgon      |HPhed      |HPtop      |Hplan      |Hrect      |Itrvl      |Line       |LineS      |Singl      |Universe   |VPgon      |VPtop      |ZeroS      |Zonot      | CP        | CPA       | CH        | CHA       |EMap       | EPM       |Itsct      |ItscA      |LiMap      | MS        | MSA       | CMS       | ReMap     | SIH       | UnionSet  | UnionSArr | Complem   |
|-------------------------------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| **Interfaces**                |           |           |           |           |           |           |           |           |           |           |           |           |           |           |   ⊎       |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |
| `LazySet`                     |   ⊎       |⊆  ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆  ⊎i      |⊆i ⊎       |⊆i ⊎  ∩    |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆  ⊎i      |⊆i ⊎i      |⊆  ⊎       |⊆i ⊎       |⊆i ⊎       |⊆i ⊎i      |⊆i ⊎i      |   ⊎       |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆  ⊎  ∩    |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎  ∩    |   ⊎  ∩    |⊆  ⊎       |
| `APolytope`                   |⊆  ⊎i      |⊆i ⊎i ∩    |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩    |⊆i ⊎i ∩    |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆  ⊎  ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `ACentrallySymmetric`         |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `ACentrallySymmetricPolytope` |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `APolygon`                    |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |-          |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `AHyperrectangle`             |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎  ∩    |⊆i ⊎i ∩i   |⊆i ⊎  ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `AHPolygon`                   |⊆i ⊎       |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩    |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |-          |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i C |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `ASingleton`                  |⊆  ⊎  ∩    |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎  ∩i   |⊆i ⊎i ∩i   |⊆  ⊎  ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎  ∩    |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i  ∩i  |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
|                               |           |           |           |           |           |           |           |           |           |           |           |           |           |           |   ⊎       |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |
| **Basic set types**           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |   ⊎       |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |
| `Ball1`                       |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Ball2`                       |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆  ⊎i ∩i   |⊆i ⊎i      |⊆ ⊎        |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆  ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆  ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `BallInf`                     |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Ballp`                       |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆  ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆  ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆  ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Ellipsoid`                   |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `EmptySet`                    |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆  ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `HalfSpace`                   |   ⊎       |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎       |⊆i ⊎i      |⊆i ⊎i ∩    |⊆i ⊎i ∩    |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎  ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `HPolygon`/`HPolygonOpt`      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i C |⊆i ⊎i ∩i C |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i C |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i ∩i C |⊆i ⊎i      |⊆i ⊎i ∩i Ci|-          |⊆i ⊎i      |⊆i ⊎i ∩i C |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i Ci|   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `HPolyhedron`                 |   ⊎       |⊆i ⊎i ∩  C |   ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|   ⊎i      |⊆i ⊎i ∩i Ci|   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i ∩    |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩    |⊆i ⊎i ∩  Ci|   ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|   ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩  Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i ∩i Ci|   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `HPolytope`                   |⊆i ⊎       |⊆i ⊎i ∩  C |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i C |⊆i ⊎i ∩i C |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩    |⊆i ⊎i ∩i C |⊆i ⊎i ∩    |⊆i ⊎i ∩  Ci|⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩  C |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i Ci|   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Hyperplane`                  |   ⊎       |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎  ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎       |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Hyperrectangle`              |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i C |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Interval`                    |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |-          |⊆i ⊎i ∩i   |-          |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |-          |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |-          |-          |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |-          |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Line`                        |   ⊎       |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |-          |   ⊎i ∩    |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎  ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `LineSegment`                 |⊆  ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |-          |⊆i ⊎i      |⊆i ⊎  ∩i   |⊆i ⊎i ∩i   |⊆  ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Singleton`                   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Universe`                    |⊆  ⊎  ∩    |⊆  ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎  ∩    |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎i ∩i   |⊆i ⊎  ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎  ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎  ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆  ⊎  ∩    |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `VPolygon`                    |⊆i ⊎i      |⊆i ⊎i ∩i C |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i C |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i ∩i Ci|-          |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i Ci|   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `VPolytope`                   |⊆i ⊎i      |⊆i ⊎i ∩i C |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩    |⊆i ⊎i ∩  C |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩  Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i Ci|   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `ZeroSet`                     |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Zonotope`                    |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i C |⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎       |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i Ci|⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i Ci|⊆i ⊎  ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i Ci|   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
|                               |           |           |           |           |           |           |           |           |           |           |           |           |           |           |   ⊎       |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |
| **Lazy set operation types**  |           |           |           |           |           |           |           |           |           |           |           |           |           |           |   ⊎       |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |           |
| `CartesianProduct`            |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `CartesianProductArray`       |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆  ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `ConvexHull`                  |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `ConvexHullArray`             |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `ExponentialMap`              |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `ExponentialProjectionMap`    |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `Intersection`                |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `IntersectionArray`           |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `LinearMap`                   |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `MinkowskiSum`                |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `MinkowskiSumArray`           |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `CacheMinkowskiSum`           |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `ResetMap`                    |   ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i      |   ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |⊆i ⊎i      |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `SymmetricIntervalHull`       |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i      |⊆i ⊎i ∩i   |   ⊎i ∩i   |   ⊎i ∩i   |⊆i ⊎i      |
| `UnionSet`                    |⊆  ⊎  ∩    |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |   ⊎       |   ⊎       |           |
| `UnionSetArray`               |⊆  ⊎  ∩    |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |⊆i ⊎i ∩i   |   ⊎       |   ⊎       |           |
| `Complement`                  |   ⊎       |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |   ⊎i      |           |           |           |

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
