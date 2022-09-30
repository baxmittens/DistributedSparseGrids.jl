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


include("../src/DistributedSparseGrids.jl")

#using UnicodePlots
#using DistributedSparseGrids
using StaticArrays 


N=5
CT = Float64
#RT = Matrix{Float64}
RT = Float64
CPType = CollocationPoint{N,CT}
HCPType = HierarchicalCollocationPoint{N,CPType,RT}

#@time begin
maxlvl = 20
nrefsteps = 8
tol = 1e-5
pointprobs = SVector{N,Int}([1 for i = 1:N])
asg = init(AHSG{N,HCPType},pointprobs)
_cpts = Set{HierarchicalCollocationPoint{N,CPType,RT}}(collect(asg))
for i = 1:nrefsteps; union!(_cpts,generate_next_level!(asg)); end
f(x::SVector{N,CT},ID::String) = 1/(sum(x.^2) + 0.3)

@time init_weights!(asg, f)

integrate(asg)
numpoints(asg)
#end

in_it = InterpolationIterator(asg,x,10)

using Distributed
addprocs(2)
works = workers()
@everywhere begin
    fun2(x,ID) = 1.0
    using StaticArrays
end 
DistributedSparseGrids.distributed_init_weights!(asg, fun2, works)
DistributedSparseGrids.integrate(asg)
using Distributed

RT = Matrix{Float64}
CPType = DistributedSparseGrids.CollocationPoint{N,CT}
HCPType = DistributedSparseGrids.HierarchicalCollocationPoint{N,CPType,RT}
asg = DistributedSparseGrids.init(DistributedSparseGrids.AHSG{N,HCPType},pointprobs)
_cpts = Set{DistributedSparseGrids.DistributedSparseGrids.HierarchicalCollocationPoint{N,CPType,RT}}(collect(asg))
for i = 1:nrefsteps; union!(_cpts,DistributedSparseGrids.generate_next_level!(asg)); end
fun3(x,ID) =ones(2,2)
DistributedSparseGrids.init_weights_inplace_ops!(asg, fun3)
DistributedSparseGrids.integrate(asg)


@everywhere fun4(x,ID) =ones(2,2)
DistributedSparseGrids.distributed_init_weights_inplace_ops!(asg, fun4, works)
DistributedSparseGrids.integrate(asg)

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