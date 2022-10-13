include("../src/DistributedSparseGrids.jl")
using Main.DistributedSparseGrids
using Test
using Distributed
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

@testset "DistributedSparseGrids.jl" begin
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
for i = 1:12
  # call generate_next_level! with tol=1e-5 and maxlevels=20
  cpts = generate_next_level!(asg, 1e-4, 12)
  # calculate weights on all worker
  distributed_init_weights!(asg, collect(cpts), fun1, ar_worker)
end

@test  isapprox(2.33905, integrate(asg), atol=1e-3)

end
