import StaticArrays: SVector, SizedArray

include("./one_dimensional_collocation_points.jl")

mutable struct CollocationPoint{N,CT<:Real,Nv}
	level::Int
	i_multi::SVector{N,Int}
	pt_idx::SVector{N,Int}
	coords::SVector{N,CT}
	interv::SVector{N,SVector{2,CT}}
	children::SVector{N,SVector{2,CollocationPoint{N,CT,Nv}}}
	parents::SVector{N,CollocationPoint{N,CT,Nv}}
	wavelet_coeffs::SVector{N,SVector{Nv,CT}}
	wavelet_cpts::SVector{N,SVector{Nv,CollocationPoint{N,CT,Nv}}}
	function CollocationPoint{N,CT,Nv}(level::Int,i_multi::SVector{N,Int},pt_idx::SVector{N,Int},coords::SVector{N,CT},interv::SVector{N,SVector{2,CT}}) where {N,CT<:Real,Nv}
		#incomplete constructor for invalid element
		cpt = new{N,CT,Nv}()
		cpt.level = -1
		cpt.i_multi = i_multi
		cpt.pt_idx = pt_idx
		cpt.coords = coords
		cpt.interv = interv
		return cpt
	end
	function CollocationPoint{N,CT,Nv}(level::Int,i_multi::SVector{N,Int},pt_idx::SVector{N,Int},coords::SVector{N,CT},interv::SVector{N,SVector{2,CT}},children::SVector{N,SVector{2,CollocationPoint{N,CT,Nv}}},parents::SVector{N,CollocationPoint{N,CT,Nv}}) where {N,CT<:Real,Nv}
		return new{N,CT,Nv}(level,i_multi,pt_idx,coords,interv,children,parents)
	end
end

include("./CollocationPoints/tuple_generation.jl")

function CollocationPoint(::Type{Val{DIM}}, parent::CollocationPoint{N,CT,Nv}, cp::CT, cp_interv::Tuple{CT,CT}, ptidx::Int, lvl_offset::Int=0) where {DIM,N,CT<:Real,Nv}
	return CollocationPoint{N,CT,Nv}(
		parent.level+lvl_offset+1,
		_unroll_(parent.i_multi,parent.i_multi[DIM]+lvl_offset+1,Val{DIM}),
		_unroll_(parent.pt_idx,ptidx,Val{DIM}),
		_unroll_(parent.coords,cp,Val{DIM}),
		_unroll_(parent.interv,SVector{2,CT}(cp_interv),Val{DIM}),
		_unroll_(SVector{N,SVector{2,CollocationPoint{N,CT,Nv}}},SVector{2,CollocationPoint{N,CT,Nv}}(invalid_collocation_point(CollocationPoint{N,CT,Nv}),invalid_collocation_point(CollocationPoint{N,CT,Nv}))),
		_unroll_(SVector{N,CollocationPoint{N,CT,Nv}},invalid_collocation_point(CollocationPoint{N,CT,Nv}),parent,Val{DIM})
	)
end

isrefined(cpt::CollocationPoint) = all(map(hasvalid,cpt.children))
getkey(cpt::CollocationPoint) = (cpt.pt_idx,cpt.i_multi)
level(cpt::CollocationPoint) = cpt.level
level(cpt::CollocationPoint, dim::Int) = cpt.i_multi[dim]
isroot(cpt::CollocationPoint{N,CT},dim::Int) where {N,CT<:Real} = cpt.pt_idx[dim] == cpt.i_multi[dim] == 1
coords(cpt::CollocationPoint) = cpt.coords
coord(cpt::CollocationPoint, dim::Int) = cpt.coords[dim]
interval(cpt::CollocationPoint, dim::Int) = cpt.interv[dim]
parent(cpt::CollocationPoint, dim::Int) = isvalid(cpt.parents[dim]) ? cpt.parents[dim] : nothing
parents(cpt::CollocationPoint) = cpt.parents
children(cpt::CollocationPoint) = isrefined(cpt) ? cpt.children : nothing
children(cpt::CollocationPoint, dim::Int) = cpt.children[dim]
isleaf(cpt::CollocationPoint) = all(map(x->!isvalid(x),children(cpt)))
intervals(cpt::CollocationPoint) = cpt.interv
Base.length(cpt::CollocationPoint) = 1
Base.isvalid(cpt::CollocationPoint) = level(cpt) > 0
hasvalid(cpts::SVector{2,CollocationPoint{N,CT,Nv}}) where {N,CT<:Real,Nv} = level(cpts[1]) > 0 || level(cpts[2]) > 0
Base.isvalid(::Nothing) = false

function Base.string(cpt::CollocationPoint{N,CT}) where {N,CT<:Real}
	str = "CollocationPoint{$N,$CT}($(cpt.level), $(cpt.coords), $(cpt.interv))"
	return str
end

function Base.display(cpt::CollocationPoint)
	print(cpt)
end

function Base.print(io::IO, cpt::CollocationPoint)
	return print(io,Base.string(cpt))
end

function Base.show(io::IO, cpt::CollocationPoint)
	return Base.print(io,cpt)
end

function root_point(::Type{CollocationPoint{N,CT,Nv}}) where {N,CT<:Real,Nv}
	rp,rp_interv = root_point(CT)
	return CollocationPoint{N,CT,Nv}(
		1,
		_unroll_(SVector{N,Int},1),
		_unroll_(SVector{N,Int},1),
		_unroll_(SVector{N,CT},rp),
		_unroll_(SVector{N,SVector{2,CT}},SVector{2,CT}(rp_interv)),
		_unroll_(SVector{N,SVector{2,CollocationPoint{N,CT,Nv}}},SVector{2,CollocationPoint{N,CT,Nv}}(invalid_collocation_point(CollocationPoint{N,CT,Nv}),invalid_collocation_point(CollocationPoint{N,CT,Nv}))),
		_unroll_(SVector{N,CollocationPoint{N,CT,Nv}},invalid_collocation_point(CollocationPoint{N,CT,Nv}))
		)
end

@generated function invalid_collocation_point(::Type{CollocationPoint{N,CT,Nv}}) where {N,CT<:Real,Nv}
	sero = zero(CT)
	cpt = CollocationPoint{N,CT,Nv}(
		-1,
		_unroll_(SVector{N,Int},-1),
		_unroll_(SVector{N,Int},-1),
		_unroll_(SVector{N,CT},sero),
		_unroll_(SVector{N,SVector{2,CT}},SVector{2,CT}(sero,sero))
		)
	return quote
		$cpt
	end
end

const CptID{N} = Tuple{SVector{N,Int},SVector{N,Int}}
const PointDict{N,CT,Nv} = Dict{CptID{N},CollocationPoint{N,CT,Nv}}

level(id::CptID{N}) where {N} = sum(id[2])-N+1

mutable struct AdaptiveSparseGrid{N,CT,Nv}
	cpts::Vector{PointDict{N,CT,Nv}}
	pointSetProperties::SVector{N,Int}
	wavelet_maxp::Int
	function AdaptiveSparseGrid{N,CT,Nv}(pointSetProperties::SVector{N,Int},wavelet_maxp::Int) where {N,CT<:Real,Nv}
		return new{N,CT,Nv}(Vector{Vector{CollocationPoint{N,CT}}}(),pointSetProperties,wavelet_maxp)
	end
	function AdaptiveSparseGrid{N,CT,Nv}(wavelet_maxp::Int) where {N,CT<:Real,Nv}
		pointSetProperties = SVector{N,Int}([1 for i = 1:N])
		return new{N,CT,Nv}(Vector{Vector{CollocationPoint{N,CT,Nv}}}(),pointSetProperties,wavelet_maxp)
	end
end
const ASG{N,CT,Nv} = AdaptiveSparseGrid{N,CT,Nv}
get_root(asg::ASG{N,CT,Nv}) where {N,CT<:Real,Nv} = first(values(asg.cpts[1]))
numlevels(asg::ASG{N,CT,Nv}) where {N,CT<:Real,Nv} = length(asg.cpts)
function Base.getindex(asg::ASG{N,CT,Nv}, ind::CptID{N}) where {N,CT,Nv}
	lvl = level(ind)
	return asg.cpts[lvl][ind]
end

function updated_interval(interv::SVector{N,SVector{2,CT}}, dim_interval::SVector{2,CT}, dim::Int) where {N,CT<:Real}
	return SVector{N,SVector{2,CT}}(interv[1:dim-1]...,dim_interval,interv[dim+1:N]...)
end

function Base.push!(asg::ASG{N,CT,Nv},cp::CollocationPoint{N,CT,Nv}) where {N,CT<:Real,Nv}
	#@assert numlevels(asg) >= (cp.level-1) && isvalid(cp)
	if cp.level > numlevels(asg)
		push!(asg.cpts,PointDict{N,CT,Nv}())
	end
	ptdict = asg.cpts[level(cp)]
	_key = getkey(cp)
	if haskey(ptdict,_key)
		ret = ptdict[_key]
		pind = -1
		for i = 1:N; if isvalid(cp.parents[i]); pind = i; break; end; end
		#parentdim = findall(map(isvalid,cp.parents))
		#@assert count(isvalid,cp.parents) == 1 && length(parentdim) == 1 && !isvalid(ret.parents[parentdim[1]]) "$(cp.parents) \n $(ret.parents) \n $(count(isvalid,cp.parents)) \n $(length(parentdim)) \n  $(isvalid(ret.parents[parentdim[1]])) \n $parentdim"
		#@assert ret.interv == updated_interval(ret.interv, cp.interv[parentdim[1]], parentdim[1])
		#ret.parents[pind] = cp.parents[pind]
		ret.parents = _unroll_dyn_dim_(ret.parents,cp.parents[pind],pind)
	else
		ret = cp
		ptdict[_key] = cp
	end
	return ret
end

function isrefined(asg::ASG, lvl::Int)
	return all(map(isrefined,values(asg.cpts[lvl])))
end

function spawn(::Type{ASG{N,CT,Nv}},pointSetProperties::SVector{N,Int},wavelet_maxp::Int) where {N,CT<:Real,Nv}
	asg = ASG{N,CT,Nv}(pointSetProperties,wavelet_maxp)
	spawn!(asg)
	return asg
end

function spawn!(asg::ASG{N,CT,Nv}) where {N,CT<:Real,Nv}
	@assert isempty(asg.cpts)
	rcp = root_point(CollocationPoint{N,CT,Nv})
	push!(asg,rcp)
	return nothing
end

function refine!(::Type{Val{DIM}},asg::ASG{N,CT,Nv},cpt::CollocationPoint{N,CT,Nv}) where {DIM,N,CT<:Real,Nv}
	tmpch1,tmpch2 = invalid_collocation_point(CollocationPoint{N,CT,Nv}),invalid_collocation_point(CollocationPoint{N,CT,Nv})
	((cl,l_interv,llvloff),(cr,r_interv,rlvloff)) = next_points(cpt.coords[DIM],cpt.interv[DIM],asg.pointSetProperties[DIM],cpt.i_multi[DIM])
	ptidx = next_level_pt_idx(cpt.pt_idx[DIM],level(cpt,DIM),asg.pointSetProperties[DIM])
	if isvalid_interv(l_interv)
		tmpch1 = CollocationPoint(Val{DIM},cpt,cl,l_interv,ptidx-1,llvloff)
		tmpch1 = push!(asg,tmpch1)
		#tmpch1.wavelet_coeffs = wavelet_coeffs(tmpch1, asg.wavelet_maxp)
	end
	if isvalid_interv(r_interv)
		tmpch2 = CollocationPoint(Val{DIM},cpt,cr,r_interv,ptidx+1,rlvloff)
		tmpch2 = push!(asg,tmpch2)
		#tmpch2.wavelet_coeffs = wavelet_coeffs(tmpch2, asg.wavelet_maxp)
	end
	return SVector{2,CollocationPoint{N,CT,Nv}}(tmpch1,tmpch2)
end

function gen_refine_code(::Type{Val{N}}) where {N}
	str = "("
	for dim = 1:N
		str*="refine!(Val{$dim},asg,cpt)"
		dim < N ? str*="," : nothing
	end
	str *= ")"
	return Meta.parse(str)
end

@generated function _unroll_refine_(asg::ASG{N,CT,Nv}, cpt::CollocationPoint{N,CT,Nv}) where {N,CT<:Real,Nv}
	refine_code = gen_refine_code(Val{N})
	return quote
		$refine_code
	end
end

function refine!(asg::ASG{N,CT,Nv}, cpt::CollocationPoint{N,CT,Nv}) where {N,CT<:Real,Nv}
	@assert !isrefined(cpt) && isvalid(cpt)
	cpt.children = SVector{N,SVector{2,CollocationPoint{N,CT,Nv}}}(_unroll_refine_(asg,cpt))
	return nothing
end

function refine!(asg::ASG, startlevel::Int=numlevels(asg), stoplevel::Int=numlevels(asg))
	@assert startlevel <= stoplevel
	startlevel == 0 ? startlevel = 1 : nothing
	refcpts = Set{CollocationPoint}()
	for i = startlevel:stoplevel
		asglvl = asg.cpts[i]
		for cpts in values(asglvl)
			if !isrefined(cpts)
				push!(refcpts,cpts)
			end
		end
	end
	for cpts in refcpts
		refine!(asg,cpts)
	end
	return nothing
end

function generate_next_level!(asg::ASG)
	asglvl = asg.cpts[end]
	for cpts in values(asglvl)
		if !isrefined(cpts)
			refine!(asg,cpts)
		end
	end
	return nothing
end


struct AncestorIterator1d{N,CT,Nv}
	cpt::CollocationPoint{N,CT,Nv}
	dim::Int
	stoplevel::Int
	function AncestorIterator1d(cpt::CollocationPoint{N,CT,Nv},dim::Int=1,stoplevel::Int=1) where {N,CT<:Real,Nv}
		return new{N,CT,Nv}(cpt,dim,stoplevel)
	end
end

function Base.length(iter::AncestorIterator1d{N,CT,Nv}) where {N,CT<:Real,Nv}
	#iter.cpt.i_multi[iter.dim]-iter.stoplevel
	len = 0
	startlevel = iter.cpt.i_multi[iter.dim]
	stoplevel = iter.stoplevel
	dim = iter.dim
	cpt = iter.cpt
	par = parent(cpt,dim)
	while isvalid(par) && len < (startlevel-stoplevel)
		len+=1
		par = parent(cpt,dim)
	end
	return len
end
Base.eltype(iter::AncestorIterator1d{N,CT,Nv}) where {N,CT<:Real,Nv} = CollocationPoint{N,CT,Nv}

function unsafe_parent(cpt::CollocationPoint{N,CT,Nv}, dim::Int) where {N,CT,Nv}
	#return cpt.parents[dim]::CollocationPoint{N,CT}
	return cpt.parents[dim]
end

function Base.iterate(iter::AncestorIterator1d{N,CT,Nv},state::CollocationPoint{N,CT,Nv}=iter.cpt) where {N,CT<:Real,Nv}
	if state.i_multi[iter.dim] > iter.stoplevel
		par = unsafe_parent(state,iter.dim)
		if isvalid(par)
			return par,par
		else
			return nothing
		end
	else
		return nothing
	end
end

struct AncestorIterator{N,CT,Nv}
	cpt::CollocationPoint{N,CT,Nv}
	stoplevel::Int
	function AncestorIterator(cpt::CollocationPoint{N,CT,Nv},stoplevel::Int=1) where {N,CT<:Real,Nv}
		return new{N,CT,Nv}(cpt,stoplevel)
	end
end
Base.eltype(iter::AncestorIterator{N,CT,Nv}) where {N,CT<:Real,Nv} = Set{CollocationPoint{N,CT,Nv}}

function Base.iterate(iter::AncestorIterator{N,CT,Nv}) where {N,CT<:Real,Nv}
	cpt = iter.cpt
	ancs = Set{CollocationPoint{N,CT,Nv}}()
	prnts = parents(cpt)
	for prntdim in prnts
		if isvalid(prntdim)
			push!(ancs, prntdim)
		end
	end
	return ancs,ancs
end

function Base.iterate(iter::AncestorIterator{N,CT,Nv},state::Set{CollocationPoint{N,CT,Nv}}) where {N,CT<:Real,Nv}
	next_ancs = Set{CollocationPoint{N,CT,Nv}}()
	if iter.stoplevel < level(first(state))
		for cpt in state
			prnts = parents(cpt)
			for prntdim in prnts
				if isvalid(prntdim)
					push!(next_ancs,prntdim)
				end
			end
		end
		if !isempty(next_ancs)
			return next_ancs,next_ancs
		else
			return nothing
		end
	else
		return nothing
	end
end


struct InterpolationIterator{N,CT,Nv}
	root_cpt::CollocationPoint{N,CT,Nv}
	x::AbstractVector{CT}
	stoplevel::Int
	function InterpolationIterator(asg::ASG{N,CT,Nv}, x::AbstractVector{CT}, stoplevel::Int=numlevels(asg)) where {N,CT<:Real,Nv}
		return new{N,CT,Nv}(get_root(asg),x,stoplevel)
	end
	function InterpolationIterator(root::CollocationPoint{N,CT,Nv}, x::AbstractVector{CT}, stoplevel::Int) where {N,CT<:Real,Nv}
		return new{N,CT,Nv}(root,x,stoplevel)
	end
end

Base.eltype(iter::InterpolationIterator{N,CT,Nv}) where {N,CT<:Real,Nv} = Set{CollocationPoint{N,CT,Nv}}

function Base.iterate(iter::InterpolationIterator{N,CT,Nv}) where {N,CT<:Real,Nv}
	root = iter.root_cpt
	next_level_iter_cpts = Set{CollocationPoint{N,CT,Nv}}()
	push!(next_level_iter_cpts,root)
	return next_level_iter_cpts,next_level_iter_cpts
end

function Base.iterate(iter::InterpolationIterator{N,CT,Nv},state::Set{CollocationPoint{N,CT,Nv}}) where {N,CT<:Real,Nv}
	next_level_iter_cpts = Set{CollocationPoint{N,CT,Nv}}()
	if iter.stoplevel > level(first(state))
		for cpt in state
			if isrefined(cpt)
				for i = 1:N
					push!(next_level_iter_cpts,unsafe_next_interpolation_descendant(cpt,iter.x,i))
				end
			end
		end
		if !isempty(next_level_iter_cpts)
			return next_level_iter_cpts,next_level_iter_cpts
		else
			return nothing
		end
	else
		return nothing
	end
end


function unsafe_next_interpolation_descendant(cpt::CollocationPoint{N,CT,Nv},x::AbstractVector{CT},dim::Int) where {N,CT<:Real,Nv}
	if isvalid(cpt.children[dim][1]) && isvalid(cpt.children[dim][2])
		if cpt.coords[dim]-x[dim]>zero(CT)
			#left child
			return cpt.children[dim][1]
		else
			#right child
			return cpt.children[dim][2]
		end
	elseif isvalid(cpt.children[dim][1])
		return cpt.children[dim][1]
	elseif isvalid(cpt.children[dim][2])
		return cpt.children[dim][2]
	else
		error()
	end
	#if cpt.coords[dim]-x[dim]>zero(CT)
	#	#left child
	#	return cpt.children[dim][1]
	#else
	#	#right child
	#	return cpt.children[dim][2]
	#end
end


import PlotlyJS
import Colors: distinguishable_colors, RGB, N0f8

function PlotlyJS.scatter(sg::Vector{PointDict{1,CT,Nv}}, lvl_offset::Bool, labelf::F=(cp->string(cp.pt_idx)*"^"*string(cp.i_multi)); kwargs...) where {CT<:Real,Nv,F<:Function}
	colors = cols = distinguishable_colors(length(sg)+1, [RGB(1,1,1)])[2:end]
	xvals = Vector{CT}()
	yvals = Vector{CT}()
	text = Vector{String}()
	clr = Vector{RGB{N0f8}}()
	for i = 1:length(sg)
		for cp in values(sg[i])
			push!(xvals,cp.coords[1])
			push!(yvals,i)
			push!(text,labelf(cp))
			push!(clr,colors[i])
		end
	end
	p = PlotlyJS.scatter(x=xvals, y=lvl_offset ? yvals : zeros(CT,length(xvals)), text=text, marker_color=clr; kwargs... )
	return p
end

function PlotlyJS.scatter(sg::Vector{PointDict{2,CT,Nv}}, lvl_offset::Bool, labelf::F=(cp->string(cp.pt_idx)*"^"*string(cp.i_multi)); kwargs...) where {CT<:Real,Nv,F<:Function}
	colors = cols = distinguishable_colors(length(sg)+1, [RGB(1,1,1)])[2:end]
	xvals = Vector{CT}()
	yvals = Vector{CT}()
	zvals = Vector{CT}()
	text = Vector{String}()
	clr = Vector{RGB{N0f8}}()
	for i = 1:length(sg)
		for cp in values(sg[i])
			push!(xvals,cp.coords[1])
			push!(yvals,cp.coords[2])
			push!(zvals,i)
			push!(text,labelf(cp))
			push!(clr,colors[i])
		end
	end
	if lvl_offset
		p = PlotlyJS.scatter3d(x=xvals, y=yvals, z=zvals, text=text, marker_color=clr; kwargs...)
	else
		p = PlotlyJS.scatter(x=xvals, y=yvals, text=text, marker_color=clr; kwargs...)
	end
	return p
end

function PlotlyJS.scatter3d(sg::Vector{PointDict{2,CT,Nv}}, labelf::F=(cp->string(cp.pt_idx)*"^"*string(cp.i_multi)); kwargs...) where {CT<:Real,Nv,F<:Function}
	colors = cols = distinguishable_colors(length(sg)+1, [RGB(1,1,1)])[2:end]
	xvals = Vector{CT}()
	yvals = Vector{CT}()
	zvals = Vector{CT}()
	text = Vector{String}()
	clr = Vector{RGB{N0f8}}()
	for i = 1:length(sg)
		for cp in values(sg[i])
			push!(xvals,cp.coords[1])
			push!(yvals,cp.coords[2])
			push!(zvals,0.0)
			push!(text,labelf(cp))
			push!(clr,colors[i])
		end
	end
	p = PlotlyJS.scatter3d(x=xvals, y=yvals, z=zvals, text=text, marker_color=clr; kwargs...)
	return p
end

function PlotlyJS.scatter(sg::Vector{PointDict{3,CT,Nv}},lvl_offset::Bool, labelf::F=(cp->string(cp.pt_idx)*"^"*string(cp.i_multi)); kwargs...) where {CT<:Real,Nv,F<:Function}
	colors = cols = distinguishable_colors(length(sg)+1, [RGB(1,1,1)])[2:end]
	xvals = Vector{CT}()
	yvals = Vector{CT}()
	zvals = Vector{CT}()
	text = Vector{String}()
	clr = Vector{RGB{N0f8}}()
	for i = 1:length(sg)
		for cp in values(sg[i])
			push!(xvals,cp.coords[1])
			push!(yvals,cp.coords[2])
			push!(zvals,cp.coords[3])
			push!(text,labelf(cp))
			push!(clr,colors[i])
		end
	end
	p = PlotlyJS.scatter3d(x=xvals, y=yvals, z=zvals, text=text, marker_color=clr; kwargs...)
	return p
end

function PlotlyJS.scatter(asg::ASG{N,CT,Nv},lvl_offset::Bool=false,labelf::F=(cp->string(cp.pt_idx)*"^"*string(cp.i_multi))) where {N,CT<:Real,Nv,F<:Function}
	PlotlyJS.scatter(asg.cpts,lvl_offset,labelf,mode="markers+text",marker_size=7,textposition="bottom center")
end

function PlotlyJS.scatter3d(asg::ASG{2,CT,Nv},labelf::F=(cp->string(cp.pt_idx)*"^"*string(cp.i_multi))) where {N,CT<:Real,Nv,F<:Function}
	PlotlyJS.scatter3d(asg.cpts,labelf,mode="markers+text",marker_size=7,textposition="bottom center")
end

include(joinpath("support","ndgrid.jl"))
function PlotlyJS.surface(asg::ASG{2,CT,Nv}, fun::F, maxp::Int, npts = 20; kwargs...) where {CT,F<:Function,Nv}
	pts = range(-1.,stop=1.,length=npts)
	xpts, ypts = ndgrid(pts,pts)
	zz = similar(xpts)
	for i = 1:npts, j=1:npts
		zz[i,j] = interpolate(asg, fun, [xpts[i,j], ypts[i,j]], maxp)
	end
	p = PlotlyJS.surface(x=xpts,y=ypts,z=zz; kwargs...)
end

function PlotlyJS.surface(asg::ASG{1,CT,Nv}, fun::F, maxp::Int, npts = 200; kwargs...) where {CT,F<:Function,Nv}
	xpts = collect(range(-1.,stop=1.,length=npts))
	ypts = similar(xpts)
	for i = 1:npts
		ypts[i] = interpolate(asg, fun, [xpts[i]], maxp)
	end
	p = PlotlyJS.scatter(x=xpts,y=ypts; kwargs...)
end

function PlotlyJS.surface(fun::F, npts = 20; kwargs...) where {F<:Function}
	pts = range(-1.,stop=1.,length=npts)
	xpts, ypts = ndgrid(pts,pts)
	zz = similar(xpts)
	for i = 1:npts, j=1:npts
		zz[i,j] = fun([xpts[i,j], ypts[i,j]])
	end
	p = PlotlyJS.surface(x=xpts,y=ypts,z=zz; kwargs...)
end

# function lagrange_basis_polynomial(x::Float64, j::Int, o::Int, xx, shift::Bool=false)
# 	@assert 0 ≤ j-1 ≤ o
# 	j0 = j-1
# 	xj = xx[j]
# 	if shift
# 		x = x + xj
# 	end
# 	if -1. <= x <= 1.
# 		res = 1.
# 		for m = 0:(j0-1)
# 			xm = xx[m+1]
# 			res *= (x-xm)/(xj-xm)
# 		end
# 		for m = (j0+1):o
# 			xm = xx[m+1]
# 			res *= (x-xm)/(xj-xm)
# 		end
# 		return res
# 	else
# 		return 0.
# 	end
# end

include("basis1d.jl")

# function lagrange_basis_polynomial(cpt::CollocationPoint{N,CT}, dim::Int, x::Float64, maxp::Int=3, level_1_constant::Bool=true) where {N, CT<:Real}
# 	#if cpt.interv[dim][1] <= x <= cpt.interv[dim][2]
# 	if maxp > 1 || (cpt.i_multi[dim] == 1 && level_1_constant)
# 		if cpt.i_multi[dim] > maxp+1
# 			anc = get_ancestors(cpt, cpt.level-maxp+2)
# 			xp = zeros(CT, length(anc)+3)
# 			for (i,an) in enumerate(anc)
# 				xp[i] = an.coords[dim]
# 			end
# 			if !isempty(anc)
# 				interv = anc[end].interv[dim]
# 			else
# 				interv = cpt.interv[dim]
# 			end
# 			xp[end-2] = interv[1]
# 			xp[end-1] = interv[2]
# 			xp[end] = cpt.coords[dim]
#
# 		else
# 			anc = get_ancestors(cpt, 1)
# 			xp = zeros(CT, length(anc)+1)
# 			for (i,an) in enumerate(anc)
# 				xp[i] = an.coords[dim]
# 			end
# 			xp[end] = cpt.coords[dim]
# 		end
# 		return lagrange_basis_polynomial(x, length(xp), length(xp)-1, xp, false)
# 	else
# 		if cpt.interv[dim][1] <= x <= cpt.interv[dim][2]
# 			if x <= cpt.coords[dim] && (cpt.coords[dim] - cpt.interv[dim][1]) > eps()
# 				res = (x - cpt.interv[dim][1]) / (cpt.coords[dim] - cpt.interv[dim][1])
# 				@assert !isnan(res) "$(cpt.coords[dim] - cpt.interv[dim][1])"
# 				return res
# 			else
# 				res = one(CT) - (x - cpt.coords[dim]) / (cpt.interv[dim][2] - cpt.coords[dim])
# 				@assert !isnan(res)
# 				return res
# 			end
# 		else
# 			return zero(CT)
# 		end
# 	end
# end

function basis_fun(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::Float64, maxp::Int=3, level_1_constant::Bool=true) where {N,CT,Nv}
	if cpt.interv[dim][1] <= x <= cpt.interv[dim][2]
		#if cpt.level > maxp
		#	return zero(CT)
		#else
		#@assert norm(lagrange_basis_polynomial(cpt, dim, x, maxp, level_1_constant) - polybasis(cpt, dim, x, maxp, level_1_constant)) <= 100*eps()
			return polybasis(cpt, dim, x, maxp, level_1_constant)
		#end
	else
		return zero(CT)
	end
end

function basis_fun(cpt::CollocationPoint{N,CT,Nv}, x::AbstractVector{CT}, maxp::Int=3, level_1_constant::Bool=true) where {N,CT<:Real,Nv}
	res = basis_fun(cpt, 1, x[1], maxp, level_1_constant)
	for dim = 2:N
		res *= basis_fun(cpt, dim, x[dim], maxp, level_1_constant)
	end
	return res
end



function wavelet_nvmoments(cpt::CollocationPoint{N,CT,Nv},dim::Int) where {N,CT,Nv}
	lvl = level(cpt,dim)
	nvmoments = min(lvl-1, Nv)
	untillvl = lvl-nvmoments
	anc = collect(AncestorIterator1d(cpt,dim,untillvl))
	return min(length(anc), nvmoments)
end

include("wavelets.jl")

function weight(asg::ASG, cpt::CollocationPoint, fun::F, maxp::Int) where {F<:Function}
	x = coords(cpt)
	fx = fun(x)
	if level(cpt) == 1
		return fx
	else
		return fx - interp_below(asg, cpt, fun, x, maxp)
	end
end

function interp_below(asg::ASG, cpt::CollocationPoint, fun::F, x::AbstractVector{T}, maxp::Int) where {T<:Real, F<:Function}
	return interpolate(asg,fun,x,maxp,level(cpt)-1)
end

function interpolate(asg::ASG{N,CT}, fun::F, x::AbstractVector{CT}, maxp::Int, stoplevel::Int=numlevels(asg)) where {N, CT, F<:Function}
	res = zero(CT) #type has to be return type of fun
	in_it = InterpolationIterator(asg,x,stoplevel)
	for cpt_set in in_it
		for cpt in cpt_set
			restmp = weight(asg, cpt, fun, maxp) * basis_fun(cpt, 1, x[1], maxp)
			for i = 2:N
				restmp *= basis_fun(cpt, i, x[i], maxp)
			end
			res += restmp
		end
	end
	return res
end

const SupportBox{N,CT} = SVector{N,SVector{2, CT}}

function intersect_support(a::T, b::T) where {CT<:Number, T<:SVector{2, CT}}
	return SVector{2,CT}((max(a[1], b[1]), min(a[2], b[2])))
end

function intersect_support(a::T, b::T) where {N, CT<:Number, T<:SupportBox{N,CT}}
	return map(intersect_support, a, b)
end

vol(sb::SupportBox{N,CT}) where {N,CT} = mapreduce(x->abs(x[2]-x[1]), *, sb)

function splitperdim(sb::SupportBox{N,CT}, dim::Int, ind::Int, which::Int) where {N,CT}
	@assert ind == 1 || ind == 2
	if which == 0
		if ind == 1
			start = sb[dim][1]
			stop = (sb[dim][2] - start) / 2 + start
		else
			intvstart = sb[dim][1]
			stop = sb[dim][2]
			start = (stop - intvstart) / 2 + intvstart
		end
	elseif which == 1
		if ind == 1
			start = sb[dim][1]
			stop = sb[dim][2]
		else
			start = sb[dim][2]
			stop = sb[dim][2]
		end
	elseif which == 2
		if ind == 1
			start = sb[dim][1]
			stop = sb[dim][1]
		else
			start = sb[dim][1]
			stop = sb[dim][2]
		end
	end
	return SVector{2,CT}((start,stop))
end

@generated function split(sb::SupportBox{N,CT}, which::SVector{N,Int}) where {N,CT}
	exprs = Vector{Expr}(undef, N+1)

	for i = 1:N
		exprs[i] = :($(Symbol("splt_$i")) = splitperdim(sb, $i, J[$i], which[$i]))
	end

	exprs[end] = Meta.parse("$(Symbol("tpl")) = ("*reduce(*,("$(Symbol("splt_$i")), " for i = 1:N))*")")
	code = Expr(:block, exprs...)
	res = quote
		dims = $(ntuple(i->2,Val{N}()))
		inds = CartesianIndices(dims)
		res = Array{SupportBox{N,CT},N}(undef, dims...)
		for J in inds
			$code
			res[J] = SupportBox{N,CT}(tpl)
		end
		return res
	end
	return res
end



function numpoints(asg::ASG{N,CT,Nv}) where {N,CT,Nv}
	res = 0
	for level in asg.cpts
		res += length(level)
	end
	return res
end

function Base.collect(asg::ASG{N,CT,Nv}) where {N,CT,Nv}
	cpts = Vector{CollocationPoint{N,CT,Nv}}(numpoints(asg))
	ctr = 0
	for level in asg.cpts
		for cpt in values(level)
			ctr += 1
			cpts[ctr] = cpt
		end
	end
	return cpts
end


mutable struct IntegrationSupport{N,M,CT}
	box::SupportBox{N,CT}
	parent::Union{Nothing,IntegrationSupport{N,M,CT}}
	children::Union{Nothing,Array{IntegrationSupport{N,CT},N}}
end

function IntegrationSupport(box::SupportBox)
	return new(box, nothing, nothing)
end



function integration_supports(asg::ASG{N,CT,Nv}) where {N, CT, Nv}
	res = Dict{Tuple{CptID{N},Int},IntegrationSupport{N,CT}}()
	rt = first(values(asg.cpts[1]))
	res[(getkey(rt),1)] = IntegrationSupport(intervals(rt))
	integration_supports!(res, rt)
	return res
end

function howtosplit(cpt::CollocationPoint{N,CT,Nv}) where {N,CT,Nv}
	x = coords(cpt)
	i = intervals(cpt)
	res = map(x, i) do y, j
		if isapprox(y, j[1], atol=10*eps())
			return 2
		elseif isapprox(y, j[2], atol=10*eps())
			return 1
		else
			return 0
		end
	end
	return res
end

function assign_parent!(child::IntegrationSupport{N,CT}, parent::IntegrationSupport{N,CT}) where {N,CT}
	child.parent = parent
	return nothing
end

function assign_children!(parent::IntegrationSupport{N,CT}, children::Array{IntegrationSupport{N,CT},N}) where {N,CT}
	parent.children = children
	return nothing
end

function integration_supports!(cpt::CollocationPoint{N,CT,Nv}, is::IntegrationSupport{N,CT}) where {N,CT,Nv}
	if isrefined(cpt)
		splt = map(IntegrationSupport, split(intervals(cpt), howtosplit(cpt)))
		assign_children!(is)
		foreach(assign_parent!, splt)
	end
	# if isrefined(cpt)
	# 	ctr = 0
	# 	for iv in split(intervals(cpt), howtosplit(cpt))
	# 		ctr += 1
	# 		if vol(iv) > 10*eps()
	# 			is = IntegrationSupport(iv)
	# 			res[(getkey(cpt),ctr)] = is
	# 			is.parent =
	# 		end
	# 	end
	for i = 1:N
		for c in children(cpt, i)
			if c !== nothing
					#integration_supports!(c, )
			end
		end
	end
	# end

	return nothing
end

function dofmap(asg::ASG{N,CT,Nv}) where {N,CT,Nv}
	ctr = 0
	npts = numpoints(asg)
	cpt_to_int = Dict{CptID{N}, Int}()
	sizehint!(cpt_to_int, npts)
	int_to_cpt = Vector{CptID{N}}(undef, npts)
	for lvl in asg.cpts
		for cptk in keys(lvl)
			ctr += 1
			cpt_to_int[cptk] = ctr
			int_to_cpt[ctr] = cptk
		end
	end
	return (cpt_to_int, int_to_cpt)
end











N = 1
#1->CPS,2->OPS,3->LOPS,4->ROPS
pointSetProperties=SVector{N,Int}([1 for i = 1:N])

maxp=3
@time begin
	asg = spawn(ASG{N,Float64,2},pointSetProperties,maxp)
	pts = collect(values(asg.cpts[end]))
	refine!(asg,pts[1])
	pts = collect(values(asg.cpts[end]))
	refine!(asg,pts[1])
	refine!(asg,pts[2])
	pts = collect(values(asg.cpts[end-1]))
	refine!(asg,pts[2])
	pts = collect(values(asg.cpts[end]))
	refine!(asg,pts[1])
	pts = collect(values(asg.cpts[end]))
	refine!(asg,pts[1])
	refine!(asg,pts[2])
end
#init_wavelet_coeffs!(asg,maxp)



import PlotlyJS
plots = PlotlyJS.GenericTrace{Dict{Symbol,Any}}[]


for lvl in asg.cpts[2:end]
	for cpt in values(lvl)
		#fun = x->basis_fun(cpt,x,maxp)
		#fun = x->wavelet(cpt,x,maxp,maxp-1)
		#display(PlotlyJS.plot(PlotlyJS.surface(asg,fun,maxp,50),PlotlyJS.Layout(title="$(getkey(cpt))")))
		#println(map(x->fun([x]),xx))
		#display(PlotlyJS.plot(PlotlyJS.surface(asg,fun,maxp,name="$(getkey(cpt))")))
		#push!(plots,PlotlyJS.surface(asg,fun,maxp,name="$(getkey(cpt))"))
	end
end

xx = linspace(-1,1,100)
for lvl in asg.cpts
	for cpt in values(lvl)
		poly = x->basis_fun(cpt,[x],maxp)
		push!(plots, PlotlyJS.scatter(x=xx, y=map(poly, xx)+level(cpt), name="$(getkey(cpt))") )
	end
end
#cpt_to_int, int_to_cpt, S, W = waveletS(asg)
#res_s = zeros(xx)
#res_w = zeros(xx)
#Ww = inv(S)'*W
#for (cid, i) in cpt_to_int
#	cpt = asg[cid]
#	for j = 1:length(xx)
#		res_s[j] += W[i] * basis_fun(cpt,xx[j:j],maxp)
#		res_w[j] += Ww[i] * wavelet(cpt,xx[j:j],maxp)
#	end
#end
#empty!(plots)
#push!(plots, PlotlyJS.scatter(x=xx, y=res_s, name="scaling") )
#push!(plots, PlotlyJS.scatter(x=xx, y=res_w, name="wavelet") )
push!(plots,PlotlyJS.scatter(asg,false,x->string(polyorder_v1(x,1,maxp))))

#PlotlyJS.plot(plots)

cpt = collect(values(asg.cpts[end]))[end]
#fun = x->wavelet(cpt,x,maxp,maxp-1)
#PlotlyJS.plot(PlotlyJS.surface(fun,100,name="$(getkey(cpt))"))
PlotlyJS.plot(plots)
