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
wasg = DistributedSparseGrids.init(DistributedSparseGrids.AHSG{N,HCPType},pointprobs)

#x =  @SArray rand(5)
#in_it = InterpolationIterator(asg,x,1)

_cpts = Set{DistributedSparseGrids.DistributedSparseGrids.HierarchicalCollocationPoint{N,CPType,RT}}(collect(wasg))
for i = 1:nrefsteps; union!(_cpts,DistributedSparseGrids.generate_next_level!(wasg)); end

fun(x,ID) = x[1]

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

using ProgressMeter
function DistributedSparseGrids.init_weights!(asg::SG, cpts::Set{HCP}, fun::F, worker_ids::Vector{Int}) where {N, HCP<:DistributedSparseGrids.AbstractHierarchicalCollocationPoint{N}, SG<:DistributedSparseGrids.AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
	@info "Starting $(length(asg)) simulatoin calls"
	hcpts = copy(cpts)
	while !isempty(hcpts)
		@info "$(length(hcpts)) collocation points remaining"
		@sync begin
			for pid in worker_ids
				if isempty(hcpts)
					break
				end
				hcpt = pop!(hcpts)
				ID = foldl((x,y)->x*"_"*y,map(string,DistributedSparseGrids.i_multi(hcpt)))*"_"*foldl((x,y)->x*"_"*y,map(string,DistributedSparseGrids.pt_idx(hcpt)))
				val = DistributedSparseGrids.coords(hcpt)
				@async begin
					_fval = remotecall_fetch(fun, pid, val, ID)
					DistributedSparseGrids.set_fval!(hcpt,_fval)
				end	
			end
		end
	end
	@info "Calculating weights"
	allasg = collect(cpts)
	for i = 1:DistributedSparseGrids.numlevels(asg)
		@info "Level $i"
		hcptar = filter(x->DistributedSparseGrids.level(x)==i,allasg)
		Threads.@threads for hcpt in hcptar	
			if DistributedSparseGrids.level(hcpt) > 1
				ret = deepcopy(DistributedSparseGrids.fval(hcpt))
				DistributedSparseGrids.interp_below!(ret,asg,hcpt)
				DistributedSparseGrids.mul!(ret,-1.0)
				DistributedSparseGrids.add!(ret,DistributedSparseGrids.fval(hcpt))
				DistributedSparseGrids.set_scaling_weight!(hcpt,ret)
			else
				DistributedSparseGrids.set_scaling_weight!(hcpt,DistributedSparseGrids.fval(hcpt))
			end
		end
	end
	return nothing
end

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