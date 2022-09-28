module DistributedSparseGrids

using StaticArrays
import StaticArrays: SVector
using Distributed
using ProgressMeter


include("./CollocationPoints.jl")

abstract type AbstractSparseGrid{N} end
abstract type AbstractHierarchicalSparseGrid{N,HCP<:AbstractHierarchicalCollocationPoint} <: AbstractSparseGrid{N} end

const PointDict{ N, HCP <: AbstractHierarchicalCollocationPoint{N}} = Dict{SVector{N,Int},Dict{SVector{N,Int},HCP}}

struct AdaptiveHierarchicalSparseGrid{N,HCP} <: AbstractHierarchicalSparseGrid{N,HCP}
	cpts::Vector{PointDict{N,HCP}}
	pointSetProperties::SVector{N,Int}
	maxp::Int
	function AdaptiveHierarchicalSparseGrid{N,HCP}(pointSetProperties::SVector{N,Int},maxp::Int) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
		return new{N,HCP}(Vector{PointDict{N,HCP}}(),pointSetProperties,maxp)
	end
end

include("./AdaptiveSparseGrids/utils.jl")

function init(::Type{AHSG{N,HCP}}, pointSetProperties::SVector{N,Int}, maxp::Int) where {N,HCP<:AbstractHierarchicalCollocationPoint{N},F<:Function}
	asg = AHSG{N,HCP}(pointSetProperties,maxp)
	init!(asg)
	return asg
end

function init!(asg::SG) where {N,HCP<:AbstractHierarchicalCollocationPoint{N},SG<:AbstractHierarchicalSparseGrid{N,HCP},F<:Function}
	@assert isempty(asg.cpts)
	rcp = root_point(HCP)
	#ID = foldl(*,map(string,i_multi(rcp)))*"_"*foldl(*,map(string,pt_idx(rcp)))
	#f = fun((coords(rcp),ID))
	#set_scaling_weight!(rcp,f)
	push!(asg,rcp)
	return nothing
end

include("./AdaptiveSparseGrids/refinement.jl")
include("./AdaptiveSparseGrids/scaling_basis.jl")
include("./AdaptiveSparseGrids/wavelet_basis.jl")
#include("./AdaptiveSparseGrids/plotting_plots.jl")

function interpolate(asg::SG, x::VCT, stoplevel::Int=numlevels(asg)) where {N,CT,VCT<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT}, HCP<:AbstractHierarchicalCollocationPoint{N,CP}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	_maxp = maxporder(asg)
	#res = zero(CT) #type has to be return type of fun
	rcp = first(asg).scaling_weight
	#RT = typeof(rcp)
	#if RT <: Number
	#	res = zero(RT)
	#elseif RT <: AbstractArray
	#	res = zeros(Float64,size(rcp))
	#else
	#	error()
	#end
	res = zero(rcp)
	tmp = zero(rcp)
	in_it = InterpolationIterator(asg,x,stoplevel)
	for cpt_set in in_it
		for hcpt in cpt_set
			#restmp = scaling_weight(hcpt)* basis_fun(hcpt, x, _maxp)
			#res += restmp
			fill!(tmp,0.0)
			add!(tmp,scaling_weight(hcpt))
			mul!(tmp,basis_fun(hcpt, x, _maxp))
			add!(res,tmp)
		end
	end
	return res
end

function interpolate!(res,asg::SG, x::VCT, stoplevel::Int=numlevels(asg)) where {N,CT,VCT<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT}, HCP<:AbstractHierarchicalCollocationPoint{N,CP}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	_maxp = maxporder(asg)
	#res = zero(CT) #type has to be return type of fun
	#rcp = first(asg).scaling_weight
	#RT = typeof(rcp)
	#if RT <: Number
	#	res = zero(RT)
	#elseif RT <: AbstractArray
	#	res = zeros(Float64,size(rcp))
	#else
	#	error()
	#end
	#res = zero(rcp)
	fill!(res,0.0)
	tmp = zero(res)
	in_it = InterpolationIterator(asg,x,stoplevel)
	for cpt_set in in_it
		for hcpt in cpt_set
			mul!(tmp,scaling_weight(hcpt),basis_fun(hcpt, x, _maxp))
			#restmp = scaling_weight(hcpt) * basis_fun(hcpt, x, _maxp)
			#for i = 2:N
			#	restmp *= basis_fun(hcpt, i, x[i], _maxp)
			#end
			#println("interp_below cpt = ",getkey(hcpt)," scaling_weight = ",scaling_weight(hcpt)," restmp = ",restmp)
			#res += restmp
			#add!(res,interp_below(asg,hcpt))
			add!(res,tmp)
		end
	end
	return nothing
end

function interpolate_wavelet(asg::SG, x::VCT, stoplevel::Int=numlevels(asg)) where {N,CT,VCT<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT}, HCP<:WaveletCollocationPoint{N,CP}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	_maxp = maxporder(asg)
	res = zero(CT) #type has to be return type of fun
	#in_it = InterpolationIterator(asg,x,stoplevel)
	#for cpt_set in in_it
		#for hcpt in cpt_set
		for hcpt in asg
			if level(hcpt) <= stoplevel
				res += hcpt.wavelet_weight * wavelet(asg, hcpt, x, _maxp)
			end
			#for i = 2:N
			#	restmp *= wavelet(hcpt, i, x[i], _maxp)
			#end
			#res += restmp
		end
	#end
	return res
end

function interp_below(asg::SG, cpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	#if !isdefined(cpt, :interp_below)
	#	cpt.interp_below = interpolate(asg,coords(cpt),level(cpt)-1)
	#end
	#return cpt.interp_below
	return interpolate(asg,coords(cpt),level(cpt)-1)
end

function interp_below!(retval, asg::SG, cpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	interpolate!(retval,asg,coords(cpt),level(cpt)-1)
	return nothing
end

function init_weights!(asg::SG, fun::F, retval_proto) where {SG<:AbstractHierarchicalSparseGrid, F<:Function}
	allasg = collect(asg)
	for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,allasg)
		Threads.@threads for hcpt in hcptar
			ID = foldl((x,y)->x*"_"*y,map(string,i_multi(hcpt)))*"_"*foldl((x,y)->x*"_"*y,map(string,pt_idx(hcpt)))
			_fval = fun(coords(hcpt),ID)
			if level(hcpt) > 1
				ret = deepcopy(retval_proto)
				interp_below!(ret,asg,hcpt)
				_fval -= ret
			end
			set_scaling_weight!(hcpt,_fval)
		end
	end
	return nothing
end

function _init_weights!(asg::SG, fun::F) where {SG<:AbstractHierarchicalSparseGrid, F<:Function}
	allasg = collect(asg)
	for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,allasg)
		Threads.@threads for hcpt in hcptar			
			_fval = fun(coords(hcpt))
			if level(hcpt) > 1
				_fval -= interp_below(asg,hcpt)
			end
			set_scaling_weight!(hcpt,_fval)
		end
	end
	return nothing
end

function __init_weights!(asg::SG, fun::F) where {SG<:AbstractHierarchicalSparseGrid, F<:Function}
	allasg = collect(asg)
	for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,allasg)
		for hcpt in hcptar			
			_fval = fun(coords(hcpt))
			if level(hcpt) > 1
				_fval -= interp_below(asg,hcpt)
			end
			set_scaling_weight!(hcpt,_fval)
		end
	end
	return nothing
end


function init_weights!(asg::SG, fun::F) where {SG<:AbstractHierarchicalSparseGrid, F<:Function}
	allasg = collect(asg)
	for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,allasg)
		Threads.@threads for hcpt in hcptar
			ID = foldl((x,y)->x*"_"*y,map(string,i_multi(hcpt)))*"_"*foldl((x,y)->x*"_"*y,map(string,pt_idx(hcpt)))
			_fval = fun(coords(hcpt),ID)
			if level(hcpt) > 1
				_fval -= interp_below(asg,hcpt)
			end
			set_scaling_weight!(hcpt,_fval)
		end
	end
	return nothing
end

function init_weights!(asg::SG, cpts::Set{HCP}, fun::F) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
	allasg = collect(cpts)
	for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,allasg)
		_hcptar = collect(hcptar)
		#Threads.@threads for hcpt in hcptar
		Threads.@threads for i = 1:length(_hcptar)
			hcpt = _hcptar[i]
			ID = foldl((x,y)->x*"_"*y,map(string,i_multi(hcpt)))*"_"*foldl((x,y)->x*"_"*y,map(string,pt_idx(hcpt)))
			_fval = fun(coords(hcpt),ID)
			set_fval!(hcpt,_fval)
			if level(hcpt) > 1
				_fval -= interp_below(asg,hcpt)
			end
			set_scaling_weight!(hcpt,_fval)
		end
	end
	return nothing
end

function init_weights!(asg::SG, cpts::Set{HCP}, fun::F, retval_proto) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
	allasg = collect(cpts)
	for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,allasg)
		_hcptar = collect(hcptar)
		#Threads.@threads for hcpt in hcptar
		Threads.@threads for i = 1:length(_hcptar)
			hcpt = _hcptar[i]
			ID = foldl((x,y)->x*"_"*y,map(string,i_multi(hcpt)))*"_"*foldl((x,y)->x*"_"*y,map(string,pt_idx(hcpt)))
			_fval = fun(coords(hcpt),ID)
			if level(hcpt) > 1
				ret = deepcopy(retval_proto)
				interp_below!(ret,asg,hcpt)
				_fval -= ret
			end
			set_scaling_weight!(hcpt,_fval)
		end
	end
	return nothing
end

#function init_weights!(asg::SG, fun::F, worker_ids::Vector{Int}) where {SG<:AbstractHierarchicalSparseGrid, F<:Function}
#	hcpts = collect(asg)
#	while !isempty(hcpts)
#		@sync begin
#			for pid in worker_ids
#				if isempty(hcpts)
#					break
#				end
#				hcpt = popfirst!(hcpts)
#				@async begin
#					ID = foldl((x,y)->x*"_"*y,map(string,i_multi(hcpt)))*"_"*foldl((x,y)->x*"_"*y,map(string,pt_idx(hcpt)))			
#					fval = remotecall_fetch(fun, pid, (coords(hcpt),ID))
#					int = interp_below(asg,hcpt)
#					set_scaling_weight!(hcpt,fval-int)
#				end
#			end
#		end
#	end
#	return nothing
#end

function init_weights!(asg::SG, cpts::Set{HCP}, fun::F, worker_ids::Vector{Int}) where {N, HCP<:AbstractHierarchicalCollocationPoint{N}, SG<:AbstractHierarchicalSparseGrid{N,HCP}, F<:Function}
	@info "Starting $(length(asg)) simulatoin calls"
	hcpts = copy(cpts)
	while !isempty(hcpts)
		@sync begin
			for pid in worker_ids
				if isempty(hcpts)
					break
				end
				hcpt = pop!(hcpts)
				ID = foldl((x,y)->x*"_"*y,map(string,i_multi(hcpt)))*"_"*foldl((x,y)->x*"_"*y,map(string,pt_idx(hcpt)))
				val = coords(hcpt)
				@async begin
					_fval = remotecall_fetch(fun, pid, val, ID)
					set_fval!(hcpt,_fval)
				end	
			end
		end
	end
	@info "Calculating weights"
	allasg = collect(cpts)
	@showprogress for i = 1:numlevels(asg)
		hcptar = filter(x->level(x)==i,allasg)
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

function max!(asg::SG) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	maxp = maxporder(asg)
	res = zero(first(asg).scaling_weight)
	add!(res,-Inf)
	for cpt in asg
		max!(res,cpt.fval)
	end
	return res
end

function min!(asg::SG) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	maxp = maxporder(asg)
	res = zero(first(asg).scaling_weight)
	add!(res,Inf)
	for cpt in asg
		min!(res,cpt.fval)
	end
	return res
end

function integrate(asg::SG) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	println("N=$N,CP=$CP,RT=$RT")
	maxp = maxporder(asg)
	res = zero(first(asg).scaling_weight)
	for cpt in asg
		res += scaling_weight(cpt) * integral_basis_fun(cpt)
	end
	return res
end

function integrate(asg::SG,fun::Function) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	println("N=$N,CP=$CP,RT=$RT")
	maxp = maxporder(asg)
	res = zero(first(asg).scaling_weight)
	for cpt in asg
		res += scaling_weight(cpt) * integral_basis_fun(cpt) * fun(coords(cpt))
	end
	return res
end

function integrate(wasg::SG,skipdims::Vector{Int}) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,HCP}}
	@assert maximum(skipdims) <= N && length(skipdims) < N
	dims = setdiff(collect(1:N), skipdims)
	_N = N-length(dims)
	_pointprobs = SVector{_N,Int}([wasg.pointSetProperties[i] for i in skipdims])
	_CPType = CollocationPoint{_N,CT}
	_HCPType = HierarchicalCollocationPoint{_N,_CPType,RT}
	asg = init(AHSG{_N,_HCPType},_pointprobs,Maxp)
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


end # module
