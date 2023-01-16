---
title: 'DistributedSparseGrids.jl: A Julia library implementing an Adaptive Sparse Grid collocation method'
tags:
  - Julia
  - stochastics
  - sparse grids
  - high-performance computing
authors:
  - name: Maximilian Bittens
    orcid: 0000-0001-9954-294X
#   equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Robert L. Gates
#   equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: Federal Institute for Geosciences and Natural Resources (BGR)
   index: 1
 - name: Zeiss
   index: 2
date: 11 October 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Abstract

Numerical integration or interpolation of high-dimensional functions is subject to the curse of dimensionality on full tensor grids. One remedy to this problem are sparse grid approximations. The additional construction effort is often times worth spending, especially for underlying functions whose evaluation is time-consuming. In the following, a Julia implementation of a local Lagrangian adaptive hierarchical sparse grid collocation method is presented, which is suitable for memory-heavy objects generated on distributed workers.

# Statement of need

[DistributedSparseGrids.jl](https://github.com/baxmittens/DistributedSparseGrids.jl) is a Julia package for integration and interpolation of functions with generic return types. There are other approaches to sparse grid approximation written in the Julia language, as [SparseGrids.jl](https://github.com/robertdj/SparseGrids.jl), [AdaptiveSparseGrids.jl](https://github.com/jacobadenbaum/AdaptiveSparseGrids.jl), [GalerkinSparseGrids.jl](https://github.com/ABAtanasov/GalerkinSparseGrids.jl) or [Tasmanian.jl](https://github.com/floswald/Tasmanian.jl). However, there is no Julia package available at the moment which is suitable if the solution of the underlying (discretized) physical problem is time and resource consuming, requiring it to be solve on either a server or cluster enviroment and, in addition, the solution is memory-heavy as well, like a Vector, Matrix, or, for example, a complete finite element solution.

# Introduction

Sparse tensor product quadrature rules, mitigating the curse of dimensionality occurring in full tensor grid constructions, were provided first by @smolyak1963quadrature. In the last two decades collocation methods were prominent in the solution of
stochastic partial differential equation as shown in @babuvska2007stochastic and @nobile2008sparse.
@ma2009adaptive were able to once again increase efficiency of the collocation approach
by introducing an error-adaptive formulation of the method, which will serve as a basis for the
collocation method described in this project. For more information about the theory of the method implemented, see e.g. @gates2015multilevel.

# Features

In the following, some key features of the implemented approach are listed.

### Arbitrary return types

[DistributedSparseGrids.jl](https://github.com/baxmittens/DistributedSparseGrids.jl) defines a ```HierarchicalCollocationPoint{N,CP,RT}``` where ```N``` is the number of dimensions, ```CP <: AbstractCollocationPoint{N,CT<:Real}```, and ```RT``` is a generic return type. ```RT``` can be conveniently defined as the type most suitable for studying the problem at hand, such as a ```Float64```, a ````Vector{Float64}```` or a ```Matrix{Float64}```, for example.<br/> Suppose the underlying physical problem stores its data in the VTU file format [@schroeder2000visualizing]. In that case, the Julia project [VTUFileHandler.jl](https://github.com/baxmittens/VTUFileHandler.jl) [@bittens2022vtufilehandler] can be used, which implements all operators needed to load complete result files into the sparse grid.

### In-place operations

Computing the weights for the hierarchical basis as well as performing interpolation and integration relies heavily on the use of *arithmetic operators*, which allocate memory. This can be a problem, especially if the result type is memory heavy. Therefore, [DistributedSparseGrids.jl](https://github.com/baxmittens/DistributedSparseGrids.jl) defines in-place variants to all of these actions given in-place variants of the arithmetic operators are defined. For further information see the [documentation](https://baxmittens.github.io/DistributedSparseGrids.jl/dev/#In-place-operations).

### Distributed computing    

If the runtime of the function to be evaluated is long, it may be necessary to distribute the load to several workers. Julia provides this functionality *out-of-the-box* via the ```Distributed``` interface. Due to the hierarchical construction of the basis and the level-wise adaptive refinement indicator, it seems necessary to include this interface in the sparse grid for a performant application of distributed computing. [DistributedSparseGrids.jl](https://github.com/baxmittens/DistributedSparseGrids.jl) uses all workers included by the ```Distributed.addprocs``` command if the ```distributed_init_weights!``` function is used to determine the hierarchical weights.

### Additional features  

- Nested one-dimensional Clenshaw-Curtis rule
- Smolyak's sparse grid construction
- local hierarchical Lagrangian basis
- different pointsets (open, closed, halfopen)
- adaptive refinement
- multi-threaded calculation of basis coefficients with ```Threads.@threads```
- integration
- experimental: integration over $X_{\sim (i)}$ (the $X_{\sim (i)}$  notation indicates the set of all variables except $X_{i}$).

# Example

Below an example of an adaptive sampling of a function with a curved singularity in 2D is provided. In \autoref{fig:example} an illustration of the sparse grid approximation is shown.

```julia
using DistributedSparseGrids
using Distributed
using StaticArrays
import PlotlyJS

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

# Sparse Grid with 4 initial levels
pp = @SVector [1,1]
asg = sparse_grid(2, pp, 4)

# add 2 worker
ar_worker = addprocs(2)

@everywhere begin
  using StaticArrays 
  # Function with curved singularity
  fun1(x::SVector{2,Float64},ID::String) =  
    (1.0-exp(-1.0*(abs(2.0 - (x[1]-1.0)^2.0 - 
      (x[2]-1.0)^2.0) +0.01)))/(abs(2-(x[1]-1.0)^2.0-(x[2]-1.0)^2.0)+0.01)
end

# calculate weights on master
init_weights!(asg, fun1)

# adaptive refine
for i = 1:20
  # call generate_next_level! with tol=1e-5 and maxlevels=20
  cpts = generate_next_level!(asg, 1e-5, 20)
  # calculate weights on all worker
  distributed_init_weights!(asg, collect(cpts), fun1, ar_worker)
end

# plot
surfplot = PlotlyJS.surface(asg, 100)
gridplot = PlotlyJS.scatter3d(asg)
PlotlyJS.plot([surfplot, gridplot])
```

![Refined sparse grid.\label{fig:example}](https://user-images.githubusercontent.com/100423479/193813765-0b7ce7b2-639a-48d3-831d-7bd5639c9fd3.PNG){height=80% width=80%}

# References