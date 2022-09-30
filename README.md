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
-	distributed function evaluation with ```julia Distributed.remotecall_fetch```
-	multi-threaded calculation the basis scaling weights  with ```julia Threads.@threads```

## Usage


### Basic usage
```julia
DistributedSparseGrids
using StaticArrays 

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
# define tolerance
tol = 1e-5
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
	union!(_cpts,generate_next_level!(asg))
end
numpoints(asg)
#define function: input are the coordinates x::SVector{N,CT} and an unique id ID::String (e.g. "1_1_1_1")
f(x::SVector{N,CT},ID::String) = sum(x.^2)

# initialize weights
@time init_weights!(asg, f)

#integrate
integrate(asg)
#interpolate
x = rand(5)*2.0.-1.0
interpolate(asg,x)
```

```julia
```

```julia
```