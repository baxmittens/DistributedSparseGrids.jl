# DistributedSparseGrids.jl

A Julia library that implements an Adaptive Sparse Grid collocation method for integrating memory-intensive objects generated on distributed workers.


## Introduction

To mitigate the "curse of dimensionality" that occurs in the integration or interpolation of high-dimensional functions using tensor-product discretizations, sparse grids use Smolyak's quadrature rule. This is particularly useful if the evaluation of the underlying function is costly. In this library an Adaptive Sparse Grid Collocation method with a local hierarchical Lagrangian basis, first proposed by [Ma and Zabaras (2010)](https://www.sciencedirect.com/science/article/pii/S002199910900028X), is implemented. For more information about the construction of Sparse Grids, see e.g. [Gates and Bittens (2015)](https://arxiv.org/abs/1509.01462).

## Install

```julia
import Pkg
```

## Implemented features

-	local hierarchical Lagrangian basis
-	adaptive refinement
-	distributed function evaluation with ```Distributed.remotecall_fetch```
-	multi-threaded calculation of basis coefficients with ```Threads.@threads```
-	usage of arbitrary return types 
-	integration
-	experimental: integration over $X_{\sim (i)}$ (the $X_{\sim (i)}$  notation indicates the set of all variables except $X_{i}$.

## Usage



### Basic usage

## Point sets

```julia
DistributedSparseGrids
using StaticArrays 

function sparse_grid(N::Int,pointprobs,nlevel=6,CT=Float64,RT=Float64)
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
pointprops = SVector{2,Int}(1,1)
asg = sparse_grid(2,pointprops) 
numpoints(asg) # returns 145
```
![image](https://user-images.githubusercontent.com/100423479/193283422-6901ef1c-e474-4a64-a143-7988c3e9be00.png)

```julia
DistributedSparseGrids
using StaticArrays 

function scalar_sparse_grid()
	# Number of dimensions
	N=5
	# Collocation point coordinate type
	CT = Float64
	# Function return type
	RT = Float64
	# define collocation point
	CPType = CollocationPoint{N,CT}
	# define hierarchical collocation point
	HCPType = HierarchicalCollocationPoint{N,CPType,RT}
	# maximum depth of grid
	maxlvl = 10
	# define refine steps
	nrefsteps = 6
	# define point properties 
	#	1->closed point set
	# 	2->open point set
	#	3->left-open point set
	#	4->right-open point set
	pointprobs = SVector{N,Int}([1 for i = 1:N])
	# init grid
	asg = init(AHSG{N,HCPType},pointprobs)
	#set of all collocation points
	cpts = Set{HierarchicalCollocationPoint{N,CPType,RT}}(collect(asg))
	# refine nrefsteps times
	for i = 1:nrefsteps
		union!(cpts,generate_next_level!(asg))
	end
	return N,CT,RT,HCPType,asg
end

N,CT,RT,HCPType,asg = scalar_sparse_grid() 

numpoints(asg) # returns 6993
#define function: input are the coordinates x::SVector{N,CT} and an unique id ID::String (e.g. "1_1_1_1")
fun1(x::SVector{N,CT},ID::String) = sum(x.^2)

# initialize weights
@time init_weights!(asg, fun1)

#integrate
integrate(asg)
#interpolate
x = rand(5)*2.0.-1.0
interpolate(asg,x)
```

## Distributed function evaluation
```julia
N,CT,RT,HCPType,asg = scalar_sparse_grid()
numpoints(asg)

using Distributed
addprocs(2)
ar_worker = workers()
@everywhere begin
	using StaticArrays
    fun2(x::SVector{N,CT},ID::String) = 1.0
end 
# Evaluate the function on 2 workers
distributed_init_weights!(asg, fun2, ar_worker)
integrate(asg)
```


```julia
```
