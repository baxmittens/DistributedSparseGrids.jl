using UnicodePlots
using DistributedSparseGrids
using StaticArrays


N=5
CT = Float64
#RT = Matrix{Float64}
RT = Float64
CPType = DistributedSparseGrids.CollocationPoint{N,CT}
HCPType = DistributedSparseGrids.HierarchicalCollocationPoint{N,CPType,RT}

maxlvl = 20
nrefsteps = 2
tol = 1e-5
pointprobs = SVector{N,Int}([1 for i = 1:N])
asg = DistributedSparseGrids.init(DistributedSparseGrids.AHSG{N,HCPType},pointprobs)

#x =  @SArray rand(5)
#in_it = InterpolationIterator(asg,x,1)

_cpts = Set{DistributedSparseGrids.DistributedSparseGrids.HierarchicalCollocationPoint{N,CPType,RT}}(collect(asg))
for i = 1:nrefsteps; union!(_cpts,DistributedSparseGrids.generate_next_level!(asg)); end

fun(x,ID) = x[1]

DistributedSparseGrids.init_weights!(asg, fun)
integrate(asg::SG)
using Distributed

addprocs(150)
worker_ids = workers()

@everywhere begin
using UnicodePlots
using DistributedSparseGrids
using StaticArrays

#fun(x,ID) = ones(10_000,10_000)
fun(x,ID) = ones(1_000,1_000)
end

N=5
CT = Float64
RT = Matrix{Float64}
CPType = DistributedSparseGrids.CollocationPoint{N,CT}
HCPType = DistributedSparseGrids.HierarchicalCollocationPoint{N,CPType,RT}

maxlvl = 20
nrefsteps = 2
tol = 1e-5
pointprobs = SVector{N,Int}([1 for i = 1:N])
wasg = DistributedSparseGrids.init!(DistributedSparseGrids.AHSG{N,HCPType},pointprobs)
_cpts = Set{DistributedSparseGrids.DistributedSparseGrids.HierarchicalCollocationPoint{N,CPType,RT}}(collect(wasg))
for i = 1:nrefsteps; union!(_cpts,DistributedSparseGrids.generate_next_level!(wasg)); end

DistributedSparseGrids.init_weights!(wasg, _cpts, fun, worker_ids)


#DistributedSparseGrids.init_weights!(wasg, _cpts, fun)

#for i = 1:nrefsteps
	#@info "refstep $i"
	#@time nchilds = generate_next_level!(wasg,tol,maxlvl)
	#@info "$(length(nchilds)) new CollocationPoints created"
	#if !isempty(nchilds)
		#@info "start init weights"
		#@time init_weights!(wasg, nchilds, fun)  
	#end
#end

#__init_weights!(wasg, fun)


# normal and inplace variant
# normal and distributed
# normal and @threads variant