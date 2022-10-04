# DistributedSparseGrids.jl

A Julia library that implements an Adaptive Sparse Grid collocation method for integrating memory-heavy objects generated on distributed workers.


## Introduction

To mitigate the "curse of dimensionality" that occurs in the integration or interpolation of high-dimensional functions using tensor-product discretizations, sparse grids use Smolyak's quadrature rule. This is particularly useful if the evaluation of the underlying function is costly. In this library an Adaptive Sparse Grid Collocation method with a local hierarchical Lagrangian basis, first proposed by [Ma and Zabaras (2010)](https://www.sciencedirect.com/science/article/pii/S002199910900028X), is implemented. For more information about the construction of Sparse Grids, see e.g. [Gates and Bittens (2015)](https://arxiv.org/abs/1509.01462).

## Install

```julia
import Pkg
```

## Implemented features

-	Nested one-dimensional Clenshaw-Curtis rule
-	Smolyak's sparse grid construction
-	local hierarchical Lagrangian basis
-	different pointsets (open, closed, halfopen)
-	adaptive refinement
-	distributed function evaluation with ```Distributed.remotecall_fetch```
-	multi-threaded calculation of basis coefficients with ```Threads.@threads```
-	usage of arbitrary return types 
-	integration
-	experimental: integration over $X_{\sim (i)}$ (the $X_{\sim (i)}$  notation indicates the set of all variables except $X_{i}$.

## Usage

### Point sets

When using sparse grids, one can choose whether the $2d$ second-level collocation points should lay on the boundary of the domain or in the middle between the origin and the boundary. (There are other choices as well.) This results in two different sparse grids, the former with almost all points on the boundary and on the coordinate axes, the latter with all points in the interior of the domain. Since one can choose for both one-dimensional children of the root point individually, there exist a multitude of different point sets for Sparse Grids.

```julia
DistributedSparseGrids
using StaticArrays 

function sparse_grid(N::Int,pointprobs,nlevel=6,RT=Float64,CT=Float64)
	# define collocation point
	CPType = CollocationPoint{N,CT}
	# define hierarchical collocation point
	HCPType = HierarchicalCollocationPoint{N,CPType,RT}
	# init grid
	asg = init(AHSG{N,HCPType},pointprobs)
	#set of all collocation points
	cpts = Set{HierarchicalCollocationPoint{N,CPType,RT}}(collect(asg))
	# fully refine grid nlevel-1 times
	for i = 1:nlevel-1
		union!(cpts,generate_next_level!(asg))
	end
	return asg
end

# define point properties 
#	1->closed point set
# 	2->open point set
#	3->left-open point set
#	4->right-open point set

asg1 = sparse_grid(2, @SVector [1,1]) 
asg2 = sparse_grid(2, @SVector [2,2]) 
asg3 = sparse_grid(2, @SVector [1,2]) 
asg4 = sparse_grid(2, @SVector [2,1]) 
asg5 = sparse_grid(2, @SVector [3,3]) 
asg6 = sparse_grid(2, @SVector [4,4]) 
asg7 = sparse_grid(2, @SVector [3,1]) 
asg8 = sparse_grid(2, @SVector [2,3]) 
asg9 = sparse_grid(2, @SVector [4,2]) 
```

<img src="https://user-images.githubusercontent.com/100423479/193283422-6901ef1c-e474-4a64-a143-7988c3e9be00.png" width="700" height="700" />

### Integration and Interpolation

```julia
asg = sparse_grid(4, @SVector [1,1,1,1]) 

#define function: input are the coordinates x::SVector{N,CT} and an unique id ID::String (e.g. "1_1_1_1")
fun1(x::SVector{N,CT},ID::String) = sum(x.^2)

# initialize weights
@time init_weights!(asg, fun1)

# integration
integrate(asg)

# interpolation
x = rand(4)*2.0 .- 1.0
val = interpolate(asg,x)	
```

### Distributed function evaluation

```julia
asg = sparse_grid(4, @SVector [1,1,1,1]) 

# add worker and register function to all workers
using Distributed
addprocs(2)
ar_worker = workers()
@everywhere begin
    using StaticArrays
    fun2(x::SVector{4,Float64},ID::String) = 1.0
end

# Evaluate the function on 2 workers
distributed_init_weights!(asg, fun2, ar_worker)
```

### Using custom return types

For custom return type ```T``` to work, following functions have to be implemented

```julia
import Base: +,-,*,/,^,zero,zeros,one,ones,copy,deepcopy

+(a::T, b::T) 
+(a::T, b::Float64) 
*(a::T, b::Float64) 
-(a::T, b::Matrix{Float64})
-(a::T, b::Float64) 
zero(a::T) 
zeros(a::T) 
one(a::T) 
one(a::T) 
copy(a::T)
deepcopy(a::T)
```

This is already the case for many data types. Below  ```RT=Matrix{Float64}``` is used.

```julia
# sparse grid with 5 dimensions and levels
pointprop = @SVector [1,2,3,4,1]
asg = sparse_grid(5, pointprop, 6, Matrix{Float64}) 

# define function: input are the coordinates x::SVector{N,CT} and an unique id ID::String (e.g. "1_1_1_1_1_1_1_1_1_1"
# for the root poin in five dimensions)
fun3(x::SVector{N,CT},ID::String) = ones(100,100).*x[1]

# initialize weights
@time init_weights!(asg, fun3)
```
### In-place operations

There are many mathematical operations executed which allocate memory while evaluting the hierarchical interpolator. Many of these allocations can be avoided by additionally implementing the ```inplace operations``` interface for data type ```T```.

```julia
import LinearAlgebra
import LinearAlgebra: mul!

DistributedSparseGrids.add!(a::T, b::T) 
DistributedSparseGrids.add!(a::T, b::Float64) 
LinearAlgebra.mul!(a::T, b::Float64) 
LinearAlgebra.mul!(a:T, b::T, c::Float64)
```

For Matrix{Float64} this interface is already implemented.

```julia
# initialize weights
@time init_weights_inplace_ops!(asg, fun3)
```

### Distributed function evaluation and in-place operations

```julia
# initialize weights
@time distributed_init_weights_inplace_ops!(asg, fun3, ar_worker)
```
