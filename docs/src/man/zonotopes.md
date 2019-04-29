# LazySets manual for Zonotope
Mathematically, a zonotope is defined as the set.
```math
Z = \\left\\{ c + ∑_{i=1}^p ξ_i g_i,~~ ξ_i \\in [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```
where ``c \\in \\mathbb{R}^n`` is its *center* and ``\\{g_i\\}_{i=1}^p``,
``g_i \\in \\mathbb{R}^n``, is the set of *generators*.
This characterization defines a zonotope as the finite Minkowski sum of line
segments.
Zonotopes can be equivalently described as the image of a unit infinity-norm
ball in ``\\mathbb{R}^n`` by an affine transformation.
## Random Zonotope
The number of generators can be controlled with the argument `num_generators`.
and dimension By `dim`.
```@example zonotope
using LazySets

Z1=rand(Zonotope,dim=1,num_generators=2)
Z2=rand(Zonotope,dim=3,num_generators=5)
```
## Plot 2D zonotopes
```@example zonotope
using Plots;gr()
Z=rand(Zonotope,dim=2,num_generators=4)
plot(Z)
```
## set membership
`∈` is the object fuction will return true if point belongs to give zonotope.
```@example zonotope
ZP = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])

[1.0, 0.2] ∈ ZP

∈([1.0, 0.1], ZP)
```
## Centre, matrix of generators, dimentions and constraints
These are the some basic properties of the zonotope can be evaluate using simple functions.
You can use center(Z) to get center but the word 'center' is there in Plots too. So do 
```@example zonotope
Z.center
Z.generators
dim(Z)
constraints_list(Z)
ngens(Z)
```
## Containment check
`⊆` Checks whether set is subset of give zonotope. 
```@example zonotope
Z = rand(Zonotope)
Z ⊆ Z
[Singleton(vi) ⊆ Z for vi in vertices_list(Z)]
Zonotope(Z.center, Z.generators/2) ⊆ Z
BallInf(zeros(2), 0.1) ⊕ Z.center ⊆ Z
BallInf(zeros(2), 1.0) ⊕ Z.center ⊆ Z
```
## Cartesian product of two zonotopes
Cartesian product create the type `CartesianProduct{Float64,Zonotope{Float64},Zonotope{Float64}}` and now we can do set operation on the object.
```@example zonotope
Z2=rand(Zonotope,dim=1,num_generators=1)
Z1=rand(Zonotope,dim=2,num_generators=1)
CartesianProduct(Z1,Z2)
```
## Minkowski-sum of zonotopes which is again a zonotope, scale and translate zonotope
The addition of generators and centre of the zonotope is the Minkowski-sum which is also zonotope.
Scale  mutiply the generators and centre by the constant.
translate adds the given vector/array with the center of zonotope.
```@example zonotope
Z1 = rand(Zonotope)
Z2 = rand(Zonotope)
plot(Z1)
plot!(Z2)
minkowski_sum(Z1,Z2) isa Zonotope
scale(4,Z1)
translate(Z1,[1.,2.])

```
## HPolytope-representation
converts the zonotope into HPolytope where it is visualised as the collection of HalfSapces
```@example zonotope
convert(HPolytope,Z1)
```
## Linear map and Affine map
Linear map checks if Dim(z)=size(Z,2),then multiply the AbstractMatrix with the  Zonotope.
Affine map does the Translation of the along the dimensions eg. X+v or X⊕v both work. 
```@example zonotope
Z = rand(Zonotope)
Z1=rand(2, 2) * Z # linear map
Z2=rand(2, 2) * Z ⊕ rand(2) # affine map
plot(Z)
plot!(Z1)
plot!(Z2)

```
## Enclosing axis-aligned box overapproximation
Is the smallest subset to contain the 
```@example zonotope
import LazySets.Approximations.ballinf_approximation
Z1=Zonotope([0.,0.],[[1.,0],[0.,1.],[1.,1.]])
Z2=Zonotope([5.,4.],[[1.,0],[0.,1.],[1.,1.]])
D=CH(Z1,Z2)
Bapprox = interval_hull(D)
plot(Z1)
plot!(Z2)
plot!(D)
plot!(Bapprox)
```
## Construction of a zonotope using line segments
The generators of the zonotpoe can be visulised as the line segments buliding the zonotope.
```@example zonotope
Z=Zonotope([0.,0.],[[1.,0.],[0.,1.]])
plot(Z)
Z1=Zonotope([0.,0.],[[1.,0],[0.,1.],[1.,1.]]) #add line joining [-[1.,1.] ,[1.,1.]]
plot!(Z1)
L1 = LineSegment(-[1., 0.], [1.0, 0.0])
plot(L1)
L2 = LineSegment(-[0., 1.], [0., 1.])
plot!(L2)
L3 = LineSegment(-[1.,1.] ,[1.,1.])
plot!(L3)
plot!(L1 ⊕ L2 ⊕ L3, 1e-3) # Zonotope
```
## Linear map of a zonotope
We can also use fuction Linearmap.
```@example zonotope
d=rand(1,2)
d = [2. 1.; 1. 3.]
ZL=Zonotope([0.,0.],[[1.,1.],[1.,0.]])
T= linear_map(d,ZL)
plot(d)
plot!(ZL)
```
## Order reduction method 
The order of a zonotope is defined as the quotient of its number of generators and its dimension.
```@example zonotope
S=rand(Zonotope,dim=2,num_generators=10)
reduce_order(S, 4)
plot(S)
E=reduce_order(S,3)
plot!(E)
reduce_order(S,2)
```
## Overapproximation of the convex hull of two zonotopes
`CH` function gives the convex hull of two zonotopes and the `overapproximate` function overapproximate the convex hull and return zonotope.
```@example zonotope
Z=Zonotope([0.,0.],[[1.,0.],[0.,1.]])
plot(Z)
Z1=Zonotope([0.,0.],[[1.,0],[0.,1.],[1.,1.]])
plot!(Z1)
D=CH(Z1,Z2)
V=overapproximate(D,Zonotope)
plot!(V)
```
##Enclosure of polytopes with zonotopes.
```@example zonotope
import LazySets.Approximations.interval_hull
Z1=rand(Zonotope,dim=2,num_generators=3)
Z2=rand(Zonotope,dim=2,num_generators=4)
Z = CH(Z1, Z2)
Happrox = interval_hull(Z)
plot(Z1)
plot!(Z2)
plot!(Happrox, alpha=0.2)
```
## The split method
split the zonotope into two. 
```@example zonotope
S=rand(Zonotope,dim=2,num_generators=10)
M,N=split(S,1)
plot(M)
plot!(N)
M,N=split(S,3)
plot(M)
plot!(N)
M,N=split(S,5)
plot(M)
plot!(N)
```
