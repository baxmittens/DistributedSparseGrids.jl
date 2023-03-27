module DistributedSparseGrids

using StaticArrays
import StaticArrays: SVector
using Distributed
using ProgressMeter

include(joinpath(".","CollocationPoints.jl"))


"""
	AbstractSparseGrid{N}

Abstract Type
	
`N` : Dimension of hierarchical sparse grid
"""
abstract type AbstractSparseGrid{N} end

"""
	AbstractHierarchicalSparseGrid{N,HCP}

Abstract Type
	
`N` : Dimension of hierarchical sparse grid
`HCP<:AbstractHierarchicalCollocationPoint` : Collocation point type
"""
abstract type AbstractHierarchicalSparseGrid{N,HCP<:AbstractHierarchicalCollocationPoint} <: AbstractSparseGrid{N} end


"""
	PointDict{ N, HCP <: AbstractHierarchicalCollocationPoint{N}}

Typedef for `Dict{SVector{N,Int},Dict{SVector{N,Int},HCP}}`
	
`N` : Dimension of hierarchical sparse grid
`HCP<:AbstractHierarchicalCollocationPoint` : Collocation point type
"""
const PointDict{ N, HCP <: AbstractHierarchicalCollocationPoint{N}} = Dict{SVector{N,Int},Dict{SVector{N,Int},HCP}}


"""
	AdaptiveHierarchicalSparseGrid{N,HCP}

Container for hierarchical collocation points 

# Fields

`cpts::Vector{PointDict{N,HCP}}` : [`DistributedSparseGrids.PointDict`](@ref) with collocation points
`pointSetProperties::SVector{N,Int}` : Vector containing all pointset properties. 
	
Pointset properties = [psp_1,...,psp_N], 
psp_i in [1,2,3,4]. 
1=>`closed point set`, 
2=>`open point set`, 
3=>`left-open point set`, 
4=>`right-open point set`.
"""
struct AdaptiveHierarchicalSparseGrid{N,HCP} <: AbstractHierarchicalSparseGrid{N,HCP}
	cpts::Vector{PointDict{N,HCP}}
	pointSetProperties::SVector{N,Int}
	function AdaptiveHierarchicalSparseGrid{N,HCP}(pointSetProperties::SVector{N,Int}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
		return new{N,HCP}(Vector{PointDict{N,HCP}}(),pointSetProperties)
	end
end

include(joinpath(".","AdaptiveSparseGrids","utils.jl"))


function init!(asg::SG) where {N,HCP<:AbstractHierarchicalCollocationPoint{N},SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	@assert isempty(asg.cpts)
	rcp = root_point(HCP)
	push!(asg,rcp)
	return nothing
end

"""
	init(::Type{AHSG{N,HCP}}, pointSetProperties::SVector{N,Int})

Initialize the sparse grid. Returns a `N`-dimensional sparse grid where only the root point has been created.


# Constructor

- `::Type{AHSG{N,HCP}}`: Define type of [`DistributedSparseGrids.ASHG`](@ref)
- `pointSetProperties::SVector{N,Int}`: Vector containing all pointset properties.

Pointset properties = [psp_1,...,psp_N], psp_i in [1,2,3,4]. 
1=>`closed point set`, 
2=>`open point set`, 
3=>`left-open point set`, 
4=>`right-open point set.`

# Example

N = 1;
pointprobs = @SVector [1];
RT = Float64;
CT = Float64;
CPType = CollocationPoint{N,CT};
HCPType = HierarchicalCollocationPoint{N,CPType,RT};
asg = init(AHSG{N,HCPType},pointprobs)

"""	
function init(::Type{AHSG{N,HCP}}, pointSetProperties::SVector{N,Int}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	asg = AHSG{N,HCP}(pointSetProperties)
	init!(asg)
	return asg
end

include(joinpath(".","AdaptiveSparseGrids","inplace_ops.jl"))
include(joinpath(".","AdaptiveSparseGrids","refinement.jl"))
include(joinpath(".","AdaptiveSparseGrids","scaling_basis.jl"))

"""
	interpolate(asg::SG, x::VCT, stplvl::Int=numlevels(asg))

Interpolate at position `x`. 

# Arguments
- `asg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}`: initialized adaptive sparse grid
- `cpts::Dict{SVector{N,Int},Dict{SVector{N,Int},HCP}}`: Dict cointaining all collocation points


"""
function interpolate(asg::SG, x::VCT, stplvl::Int=numlevels(asg)) where {N,CT,VCT<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT}, HCP<:AbstractHierarchicalCollocationPoint{N,CP}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	rcp = scaling_weight(first(asg))
	res = zero(rcp)
	in_it = InterpolationIterator(asg,x,stplvl)
	for hcpt in in_it
		res += scaling_weight(hcpt) * basis_fun(hcpt, x, 1)
	end
	return res
end

# Recursion can arrive at Colloaction Points several times. How to avoid that?

#function interpolate_recursive(asg::SG, hcpt::HCP, x::VCT, stplvl::Int=numlevels(asg)) where {N,CT,VCT<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT}, HCP<:AbstractHierarchicalCollocationPoint{N,CP}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
#	rcp = scaling_weight(hcpt)
#	res = zero(rcp)
#	bf = basis_fun(hcpt, x, 1)
#	#if bf > 0.0
#		println(hcpt)
#		res += scaling_weight(hcpt) .* bf
#		if level(hcpt) < stplvl && isrefined(hcpt)
#			for dim = 1:N
#				ncp = next_interpolation_descendant(hcpt,x[dim],dim)	
#				res += interpolate_recursive(asg, ncp, x, stplvl)
#			end
#		end
#	#end
#	return res
#end

#function interpolate_recursive(asg::SG, x::VCT, stplvl::Int=numlevels(asg)) where {N,CT,VCT<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT}, HCP<:AbstractHierarchicalCollocationPoint{N,CP}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
#	root = get_root(asg)
#	return interpolate_recursive(asg, root, x, stplvl) 
#end
#
#function interpolate_recursive!(res::RT, tmp::RT, asg::SG, hcpt::HCP, x::VCT, stplvl::Int=numlevels(asg)) where {N,RT,CT,VCT<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT}, HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
#	bf = basis_fun(hcpt, x, 1)
#	if bf > 0.0
#		mul!(tmp,scaling_weight(hcpt),bf)
#		add!(res,tmp)
#		if level(hcpt) < stplvl #&& isrefined(hcpt)
#			for dim = 1:N
#				ncp = next_interpolation_descendant(hcpt,x[dim],dim)	
#				interpolate_recursive!(res, tmp, asg, ncp, x, stplvl)
#			end
#		end
#	end
#	return nothing
#end

#function interpolate_recursive!(res::RT, asg::SG, x::VCT, stplvl::Int=numlevels(asg)) where {N,RT,CT,VCT<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT}, HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
#	fill!(res,0.0)
#	tmp = deepcopy(res)
#	root = get_root(asg)
#	interpolate_recursive!(res, tmp, asg, root, x, stplvl) 
#	return nothing
#end

function interpolate!(res::RT, asg::SG, x::VCT, stplvl::Int=numlevels(asg)) where {N,RT,CT,VCT<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT}, HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	fill!(res,0.0)
	tmp = zero(res)
	in_it = InterpolationIterator(asg,x,stplvl)
	#for cpt_set in in_it
	for hcpt in in_it
		#for hcpt in cpt_set
		mul!(tmp,scaling_weight(hcpt),basis_fun(hcpt, x, 1))
		add!(res,tmp)
		#end
	end
	return nothing
end

function interp_below!(retval::RT, asg::SG, cpt::HCP) where {N,RT,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	interpolate!(retval,asg,coords(cpt),level(cpt)-1)
	#interpolate_recursive!(retval,asg,coords(cpt),level(cpt)-1)
	return nothing
end

function interp_below(asg::SG, cpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	return interpolate(asg,coords(cpt),level(cpt)-1)
	#return interpolate_recursive(asg,coords(cpt),level(cpt)-1)
end

"""
	init_weights!(asg::SG, cpts::AbstractVector{HCP}, fun::F)

(Re-)Computes all weights in `cpts`. 

# Arguments
- `asg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}`: adaptive sparse grid
- `cpts::AbstractVector{HCP}`: all weights of the collocation points in `cpts` will be (re-)calculated.
- `fun`::Function to be interpolated.

"""
function init_weights!(asg::SG, cpts::AbstractVector{HCP}, fun::F) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid, F<:Function}
	for i = 1:numlevels(asg)	# do it level-wise since interp_below operates on the l-1-level interpolator
		#println("level $i")
		hcptar = filter(x->level(x)==i,cpts)
		Threads.@threads for hcpt in hcptar
		#for hcpt in hcptar
			ID = idstring(hcpt)
			_fval = fun(coords(hcpt),ID)
			set_fval!(hcpt,_fval)
			if level(hcpt) > 1
				sw = _fval - interp_below(asg,hcpt)
				set_scaling_weight!(hcpt,sw)
			else
				set_scaling_weight!(hcpt,_fval)
			end
			
		end
	end
end

"""
	init_weights!(asg::SG, fun::F)

(Re-)Computes all weights in `asg`. 

# Arguments
- `asg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}`: adaptive sparse grid
- `fun`::Function to be interpolated.

"""
function init_weights!(asg::SG, fun::F) where {SG<:AbstractHierarchicalSparseGrid, F<:Function}
	allasg = collect(asg)
	init_weights!(asg, allasg, fun)
	return nothing
end

"""
	init_weights_inplace_ops!(asg::SG, cpts::AbstractVector{HCP}, fun::F)

(Re-)Computes all weights in `cpts` with in-place operations. In-place functions `mul!(::RT,::RT)`,`mul!(::RT,::Float64)`,`add!(::RT,::RT)` have to be defined.

# Arguments
- `asg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}`: adaptive sparse grid
- `cpts::AbstractVector{HCP}`: all weights of the collocation points in `cpts` will be (re-)calculated.
- `fun`::Function to be interpolated.

"""
function init_weights_inplace_ops!(asg::SG, cpts::AbstractVector{HCP}, fun::F) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid, F<:Function}
	for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,cpts)
		Threads.@threads for hcpt in hcptar
			ID = idstring(hcpt)
			_fval = fun(coords(hcpt),ID)
			set_fval!(hcpt,_fval)
			if level(hcpt) > 1
				scalweight = deepcopy(_fval)
				interp_below!(scalweight,asg,hcpt)
				mul!(scalweight,-1.0)
				add!(scalweight,fval(hcpt))
				set_scaling_weight!(hcpt,scalweight)
			else
				set_scaling_weight!(hcpt,fval(hcpt))
			end
		end
	end
end

"""
	init_weights_inplace_ops!(asg::SG, fun::F)

(Re-)Computes all weights in `asg` with in-place operations. In-place functions `mul!(::RT,::RT)`,`mul!(::RT,::Float64)`,`add!(::RT,::RT)` have to be defined.

# Arguments
- `asg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}`: adaptive sparse grid



"""
function init_weights_inplace_ops!(asg::SG, fun::F) where {SG<:AbstractHierarchicalSparseGrid, F<:Function}
	allasg = collect(asg)
	init_weights_inplace_ops!(asg, allasg, fun)
	return nothing
end

#function distributed_fvals!(asg::SG, cpts::AbstractVector{HCP}, fun::F, worker_ids::Vector{Int}) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
#@info "Starting $(length(cpts)) simulatoin calls"
#	hcpts = copy(cpts)
#	while !isempty(hcpts)
#		@sync begin
#			for pid in worker_ids
#				if isempty(hcpts)
#					break
#				end
#				hcpt = pop!(hcpts)
#				ID = idstring(hcpt)
#				val = coords(hcpt)
#				@async begin
#					_fval = remotecall_fetch(fun, pid, val, ID)
#					set_fval!(hcpt,_fval)
#				end	
#			end
#		end
#	end
#	return nothing
#end

function distributed_fvals!(asg::SG, cpts::AbstractVector{HCP}, fun::F, worker_ids::Vector{Int}) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
	@info "Starting $(length(cpts)) simulatoin calls"
	#wp = WorkerPool(worker_ids);
	wp = WorkerPool(workers())
	@sync begin
		for hcpt in cpts
			@async begin
				#println("async remotecall")
				ID = idstring(hcpt)
				val = coords(hcpt)
				_fval = remotecall_fetch(fun, wp, val, ID)
				set_fval!(hcpt,_fval)
			end
		end
	end
	return nothing
end


"""
	distributed_init_weights!(asg::SG, cpts::AbstractVector{HCP}, fun::F, worker_ids::Vector{Int})

Computes all weights in `cpts` on all workers in `worker_ids`. 

# Arguments
- `asg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}`: adaptive sparse grid
- `cpts::AbstractVector{HCP}`: all weights of the collocation points in `cpts` will be (re-)calculated.
- `fun`::Function to be interpolated.
- `worker_ids`: All available workers (can be added via `using Distributed; addprocs(...)`).

"""
function distributed_init_weights!(asg::SG, cpts::AbstractVector{HCP}, fun::F, worker_ids::Vector{Int}) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
	distributed_fvals!(asg, cpts, fun, worker_ids)
	@info "Calculating weights"
	for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,cpts)
		Threads.@threads for hcpt in hcptar
			ID = idstring(hcpt)
			if level(hcpt) > 1
				ret = fval(hcpt) - interp_below(asg,hcpt)
				set_scaling_weight!(hcpt,ret)
			else
				set_scaling_weight!(hcpt,fval(hcpt))
			end
		end
	end
	return nothing
end

"""
	distributed_init_weights!(asg::SG, fun::F, worker_ids::Vector{Int})

Computes all weights in `asg` on all workers in `worker_ids`. 

# Arguments
- `asg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}`: adaptive sparse grid
- `fun`::Function to be interpolated.
- `worker_ids`: All available workers (can be added via `using Distributed; addprocs(...)`).

"""
function distributed_init_weights!(asg::SG, fun::F, worker_ids::Vector{Int}) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
	allasg = collect(asg)
	distributed_init_weights!(asg, allasg, fun, worker_ids)
	return nothing
end

"""
	distributed_init_weights_inplace_ops!(asg::SG, cpts::AbstractVector{HCP}, fun::F, worker_ids::Vector{Int})

Computes all weights in `cpts` on all workers in `worker_ids`. In-place functions `mul!(::RT,::RT)`,`mul!(::RT,::Float64)`,`add!(::RT,::RT)` have to be defined.

# Arguments
- `asg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}`: adaptive sparse grid
- `cpts::AbstractVector{HCP}`: all weights of the collocation points in `cpts` will be (re-)calculated.
- `fun`::Function to be interpolated.
- `worker_ids`: All available workers (can be added via `using Distributed; addprocs(...)`).

"""
function distributed_init_weights_inplace_ops!(asg::SG, cpts::AbstractVector{HCP}, fun::F, worker_ids::Vector{Int}) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
	distributed_fvals!(asg, cpts, fun, worker_ids)
	@info "Calculating weights"
	@showprogress for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,cpts)
		Threads.@threads for hcpt in hcptar	
			if level(hcpt) > 1
				ret = deepcopy(fval(hcpt))
				interp_below!(ret,asg,hcpt)
				mul!(ret,-1.0)
				add!(ret,fval(hcpt))
				set_scaling_weight!(hcpt,ret)
			else
				set_scaling_weight!(hcpt,fval(hcpt))
			end
		end
	end
	return nothing
end

"""
	distributed_init_weights_inplace_ops!(asg::SG, fun::F, worker_ids::Vector{Int})

(Re-)computes all weights in `asg` on all workers in `worker_ids`. In-place functions `mul!(::RT,::RT)`,`mul!(::RT,::Float64)`,`add!(::RT,::RT)` have to be defined.

# Arguments
- `asg::SG<:AbstractHierarchicalSparseGrid{N,HCP}}`: adaptive sparse grid
- `fun`::Function to be interpolated.
- `worker_ids`: All available workers (can be added via `using Distributed; addprocs(...)`).

"""
function distributed_init_weights_inplace_ops!(asg::SG, fun::F, worker_ids::Vector{Int}) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
	allasg = collect(asg)
	distributed_init_weights_inplace_ops!(asg, allasg, fun, worker_ids)
	return nothing
end

function max!(asg::SG) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	res = zero(first(asg).scaling_weight)
	add!(res,-Inf)
	for cpt in asg
		max!(res,cpt.fval)
	end
	return res
end

function min!(asg::SG) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	res = zero(first(asg).scaling_weight)
	add!(res,Inf)
	for cpt in asg
		min!(res,cpt.fval)
	end
	return res
end

"""
    integrate(asg::SG) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}

Integrates the sparse grid. Returns an instance of RT.

# Constructor
- `asg::SG`: adaptie sparse grid

"""	
function integrate(asg::SG) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	println("N=$N,CP=$CP,RT=$RT")
	res = zero(first(asg).scaling_weight)
	for cpt in asg
		res += scaling_weight(cpt) * integral_basis_fun(cpt)
	end
	return res
end

function integrate_inplace_ops(asg::SG) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	println("N=$N,CP=$CP,RT=$RT")
	res = zero(first(asg).scaling_weight)
	tmp = zero(first(asg).scaling_weight)
	for cpt in asg
		bf = integral_basis_fun(cpt)
		sw = scaling_weight(cpt)
		mul!(tmp, sw, bf)
		add!(res, tmp) 
	end
	return res
end

#function integrate(asg::SG,fun::Function) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
#	println("N=$N,CP=$CP,RT=$RT")
#	maxp = maxporder(asg)
#	res = zero(first(asg).scaling_weight)
#	for cpt in asg
#		res += scaling_weight(cpt) * integral_basis_fun(cpt) * fun(coords(cpt))
#	end
#	return res
#end

function integrate(wasg::SG,skipdims::Vector{Int}) where {N,RT,CT,CP<:AbstractCollocationPoint{N,CT}, HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	@assert maximum(skipdims) <= N && length(skipdims) < N
	dims = setdiff(collect(1:N), skipdims)
	_N = N-length(dims)
	_pointprobs = SVector{_N,Int}([wasg.pointSetProperties[i] for i in skipdims])
	_CPType = CollocationPoint{_N,CT}
	_HCPType = HierarchicalCollocationPoint{_N,_CPType,RT}
	asg = init(AHSG{_N,_HCPType},_pointprobs)
	__cpts = Set{HierarchicalCollocationPoint{_N,_CPType,RT}}(collect(asg))
	pda = Vector{Dict{SVector{_N,Int},Dict{SVector{_N,Int},RT}}}()
	push!(pda, Dict{SVector{_N,Int},Dict{SVector{_N,Int},RT}}())
	for i = 1:length(wasg.cpts)-1
		union!(__cpts,generate_next_level!(asg))
		push!(pda, Dict{SVector{_N,Int},Dict{SVector{_N,Int},RT}}())
	end
	rpsw = scaling_weight(first(wasg))
	for hcpt in asg
		set_scaling_weight!(hcpt, zero(rpsw))
	end
	for hcpt in wasg
		res = one(CT)
		for d in dims
			res *= integral_basis_fun(hcpt, d)
		end
		res *= scaling_weight(hcpt)
		_lvl = sum(i_multi(hcpt)[skipdims])-_N+1
		if haskey(pda[_lvl],i_multi(hcpt)[skipdims]) 
			if haskey(pda[_lvl][i_multi(hcpt)[skipdims]],pt_idx(hcpt)[skipdims])
				pda[_lvl][i_multi(hcpt)[skipdims]][pt_idx(hcpt)[skipdims]] += res
			else
				pda[_lvl][i_multi(hcpt)[skipdims]][pt_idx(hcpt)[skipdims]] = res
			end
		else
			pda[_lvl][i_multi(hcpt)[skipdims]] = Dict{SVector{_N,Int},RT}()
			pda[_lvl][i_multi(hcpt)[skipdims]][pt_idx(hcpt)[skipdims]] = res
		end
	end
	for (_lvl,_dict1) in enumerate(pda)
		for (_i_multi,_dict2) in _dict1
			for (_pt_idx,sw) in _dict2
				#println(_lvl," ",_i_multi," ",_pt_idx," ",sw)
				set_scaling_weight!(asg.cpts[_lvl][_i_multi][_pt_idx],sw)
			end
		end 
	end
	return asg
end

include(joinpath(".","AdaptiveSparseGrids","plotting.jl"))
export CollocationPoint, HierarchicalCollocationPoint, AHSG, init, generate_next_level!, init_weights!, distributed_init_weights!, init_weights_inplace_ops!, distributed_init_weights_inplace_ops!, integrate, interpolate, interpolate!, integrate_inplace_ops

end # module
