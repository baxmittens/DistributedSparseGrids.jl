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

# Define normal distribution
norml(x,μ,σ) = 1/(σ*sqrt(2*pi))*exp(-0.5*((x[1]-μ)/σ)^2)

# Define 3 1d sparse grids with different points sets
asg11 = sparse_grid(1, @SVector [1]) 
asg12 = sparse_grid(1, @SVector [2]) 
asg13 = sparse_grid(1, @SVector [3]) 

# 1d Normal distribution with μ=0.0 and σ=0.1
norml1d(x::SVector{1,Float64},ID::String) = norml(x[1],0.0,0.1)

# Init sparse-grids
init_weights!(asg11, norml1d)
init_weights!(asg12, norml1d)
init_weights!(asg13, norml1d)

# adaptive refine
while true
  # call generate_next_level! with tol=1e-4 and maxlevels=30
  cpts11 = generate_next_level!(asg11, 1e-4, 30)
  cpts12 = generate_next_level!(asg12, 1e-4, 30)
  cpts13 = generate_next_level!(asg13, 1e-4, 30)
  #
  if isempty(cpts11) && isempty(cpts12) && isempty(cpts13)
    break
  end
  # calculate weights
  init_weights!(asg11, collect(cpts11), norml1d)
  init_weights!(asg12, collect(cpts12), norml1d)
  init_weights!(asg13, collect(cpts13), norml1d)
end

# integrating over pdf should always result in 1.0
@test  isapprox(1.0, integrate(asg11), atol=1e-4)
@test  isapprox(1.0, integrate(asg12), atol=1e-4)
@test  isapprox(1.0, integrate(asg13), atol=1e-4)

# Define 2 2d sparse grids with different points sets
asg21 = sparse_grid(2, @SVector [1,1]) 
asg22 = sparse_grid(2, @SVector [2,2]) 

# 2d Normal distribution with μ=0.0 and σ=0.1
norml2d(x::SVector{2,Float64},ID::String) = norml(x[1],0.0,0.1)*norml(x[2],0.0,0.1)

# init weights
init_weights!(asg21, norml2d)
init_weights!(asg22, norml2d)

# adaptive refine
while true
  # call generate_next_level! with tol=1e-4 and maxlevels=30
  cpts21 = generate_next_level!(asg21, 1e-4, 30)
  cpts22 = generate_next_level!(asg22, 1e-4, 30)
  #
  if isempty(cpts21) && isempty(cpts22)
    break
  end
  # calculate weights
  init_weights!(asg21, collect(cpts21), norml2d)
  init_weights!(asg22, collect(cpts22), norml2d)
end

# integrating over pdf should always result in 1.0
@test  isapprox(1.0, integrate(asg21), atol=1e-4)
@test  isapprox(1.0, integrate(asg22), atol=1e-4)

# Define 1 3d sparse grids
asg31 = sparse_grid(3, @SVector [1,1,1]) 

# Normal distribution with μ=0.0 and σ=0.1
norml3d(x::SVector{3,Float64},ID::String) = norml(x[1],0.0,0.1)*norml(x[2],0.0,0.1)*norml(x[3],0.0,0.1)

init_weights!(asg31, norml3d)

# adaptive refine
while true
  # call generate_next_level! with tol=1e-4 and maxlevels=30
  cpts31 = generate_next_level!(asg31, 1e-4, 30)
  #
  if isempty(cpts31)
    break
  end
  # calculate weights
  init_weights!(asg31, collect(cpts31), norml3d)
end

# integrating over pdf should always result in 1.0
@test  isapprox(1.0, integrate(asg31), atol=1e-4)

# Sparse Grid with 4 initial levels
pp = @SVector [1,1]
# The quality of the error estimator depends on the number of colloaction points. 
# For very few collocation points the error estimator is not always reliable.
initial_levels = 4 
asg = sparse_grid(2, pp, initial_levels)

# add 2 worker
ar_worker = addprocs(2)

@everywhere begin
  using StaticArrays 
  # Function with curved singularity. This is an example of the most complex function the sparse grid is able to sample.
  # Function with jumps are not sampable due to the lagrangian basis
  fun1(x::SVector{2,Float64},ID::String) =  
    (1.0-exp(-1.0*(abs(2.0 - (x[1]-1.0)^2.0 - 
      (x[2]-1.0)^2.0) +0.01)))/(abs(2-(x[1]-1.0)^2.0-(x[2]-1.0)^2.0)+0.01)
end

# calculate weights on master
init_weights!(asg, fun1)

# adaptive refine
while true
  # call generate_next_level! with tol=1e-4 and maxlevels=30
  cpts = generate_next_level!(asg, 1e-4, 30)
  #
  if isempty(cpts)
    break
  end
  # calculate weights on all worker
  distributed_init_weights!(asg, collect(cpts), fun1, ar_worker)
end

analytical_integral = 2.33905
# due to wolframalpha see
# https://www.wolframalpha.com/input?i=integrate+%281.0-exp%28-1.0*%28abs%282.0+-+%28x-1.0%29%5E2.0+-+%28y-1.0%29%5E2.0%29+%2B0.01%29%29%29%2F%28abs%282-%28x-1.0%29%5E2.0-%28y-1.0%29%5E2.0%29%2B0.01%29+x%3D-1+to+1+%2C+y%3D-1+to+1

@test  isapprox(analytical_integral, integrate(asg), atol=1e-4)

end
