cpt(hcpt::AbstractHierarchicalCollocationPoint) = hcpt.cpt
i_multi(hcpt::AbstractHierarchicalCollocationPoint) = i_multi(hcpt.cpt)
i_multi(hcpt::AbstractHierarchicalCollocationPoint,dim::Int) = i_multi(hcpt.cpt,dim)
pt_idx(hcpt::AbstractHierarchicalCollocationPoint) = pt_idx(hcpt.cpt)
pt_idx(hcpt::AbstractHierarchicalCollocationPoint,dim::Int) = pt_idx(hcpt.cpt,dim)
coords(hcpt::AbstractHierarchicalCollocationPoint) = coords(hcpt.cpt)
coord(hcpt::AbstractHierarchicalCollocationPoint,dim::Int) = coord(hcpt.cpt,dim)
intervals(hcpt::AbstractHierarchicalCollocationPoint) = intervals(hcpt.cpt)
interval(hcpt::AbstractHierarchicalCollocationPoint,dim::Int) = interval(hcpt.cpt,dim)

getkey(hcpt::HCP) where {HCP<:AbstractHierarchicalCollocationPoint} = (pt_idx(hcpt),i_multi(hcpt))
level(hcpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}} = sum(i_multi(hcpt))-N+1
level(hcpt::HCP, dim::Int) where {HCP<:AbstractHierarchicalCollocationPoint} = i_multi(hcpt,dim)
isroot(hcpt::HCP,dim::Int) where {HCP<:AbstractHierarchicalCollocationPoint} = pt_idx(hcpt,dim) == i_multi(hcpt,dim) == 1

Base.length(hcpt::HCP) where {HCP<:AbstractHierarchicalCollocationPoint} = 1
Base.isvalid(hcpt::HCP) where {HCP<:AbstractHierarchicalCollocationPoint} = level(hcpt) > 0
Base.isvalid(::Nothing) = false
hasvalid(hcpts::SVector{2,HCP}) where {HCP<:AbstractHierarchicalCollocationPoint} = level(hcpts[1]) > 0 || level(hcpts[2]) > 0
Base.hash(hcpt::HCP) where {HCP<:AbstractHierarchicalCollocationPoint} = hash(getkey(hcpt))


parent(hcpt::AbstractHierarchicalCollocationPoint, dim::Int) = hcpt.parents[dim]
parents(hcpt::AbstractHierarchicalCollocationPoint) = hcpt.parents
children(hcpt::AbstractHierarchicalCollocationPoint) = hcpt.children
children(hcpt::AbstractHierarchicalCollocationPoint, dim::Int) = hcpt.children[dim]
child(hcpt::AbstractHierarchicalCollocationPoint, dim::Int, i::Int) = hcpt.children[dim][i]
isleaf(hcpt::AbstractHierarchicalCollocationPoint) = all(map(x->!isvalid(x),children(hcpt)))

isrefined(hcpt::AbstractHierarchicalCollocationPoint,dim::Int) = i_multi(hcpt,dim) == 2 ? hasvalid(children(hcpt,dim)) : all(map(isvalid,children(hcpt,dim)))
isrefined(hcpt::AbstractHierarchicalCollocationPoint{N}) where {N} = all(map(x->isrefined(hcpt,x),1:N))

setparents!(hcpt::HCP, parents::SVector{N,HCP}) where {N,HCP <: AbstractHierarchicalCollocationPoint{N}} = begin hcpt.parents = parents; return nothing; end
setchildren!(hcpt::HCP, children::SVector{N,SVector{2,HCP}}) where {N,HCP <: AbstractHierarchicalCollocationPoint{N}} = begin hcpt.children = children; return nothing; end


idstring(hcpt::HCP) where {N,HCP <: AbstractHierarchicalCollocationPoint{N}} = foldl((x,y)->x*"_"*y,map(string,i_multi(hcpt)))*"_"*foldl((x,y)->x*"_"*y,map(string,pt_idx(hcpt))) 


function gen_set_child_code(::Type{Val{N}},::Type{Val{DIM}},::Type{Val{I}}) where {N,DIM,I}
	str = "("
	for dim = 1:N
		if dim == DIM
			str*= "_unroll_(children(hcpt,$dim),chld,Val{I})"
		else
			str*="children(hcpt,$dim)"
		end
		dim < N ? str*="," : nothing
	end
	str *= ")"
	return Meta.parse(str)
end

@generated function setchild!(hcpt::HCP, ::Type{Val{DIM}}, ::Type{Val{I}}, chld::HCP) where {N,DIM,I,HCP <: AbstractHierarchicalCollocationPoint{N}}
	code = gen_set_child_code(Val{N}, Val{DIM}, Val{I})
	return quote
		hcpt.children = SVector{N,SVector{2,$HCP}}($code)
	end
end

scaling_weight(hcpt::HCP) where {HCP <: AbstractHierarchicalCollocationPoint} = hcpt.scaling_weight
set_scaling_weight!(hcpt::HCP,weight::RT) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}} = hcpt.scaling_weight = weight

set_fval!(hcpt::HCP,fval::RT) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}} = hcpt.fval = fval
fval(hcpt::HCP) where {HCP <: AbstractHierarchicalCollocationPoint} = hcpt.fval

#import Base.convert
#function convert(::Type{HCP1}, hcpt::HCP2) where {N,HCP1 <: AbstractHierarchicalCollocationPoint{N}, HCP2 <: AbstractHierarchicalCollocationPoint{N}}
#	println(HCP1)
#	println(HCP2)
#	return HCP1(cpt(hcpt),children(hcpt),parents(hcpt),scaling_weight(hcpt))
#end

function Base.string(hcpt::HCP) where {HCP <: AbstractHierarchicalCollocationPoint}
	str = "$HCP($(cpt(hcpt)))"
	return str
end

function Base.display(hcpt::AbstractHierarchicalCollocationPoint)
	print(hcpt)
end

function Base.print(io::IO, hcpt::AbstractHierarchicalCollocationPoint)
	return print(io,Base.string(hcpt))
end

function Base.show(io::IO, hcpt::AbstractHierarchicalCollocationPoint)
	return Base.print(io,hcpt)
end

#const CptID{N} = Tuple{SVector{N,Int},SVector{N,Int}}

function possible_parent_cptids(::Type{Val{DIM}},id::CptID{N},psp::Int) where {DIM,N}
	imulti = id[2]
	ptidx = id[1]		
	ptidx1d = unsafe_parent_pt_idx(ptidx[DIM],imulti[DIM],psp)
	imulti1d = imulti[DIM]-1
	_imulti = _unroll_(imulti,imulti1d,Val{DIM})
	_ptidx = _unroll_(ptidx,ptidx1d,Val{DIM})
	return (_ptidx,_imulti)
end


function gen_possible_parent_cptids_code(::Type{Val{N}}) where {N}
		str = "("
	for dim = 1:N
		str*="possible_parent_cptids(Val{$dim},id,psp[$dim])"
		dim < N ? str*="," : nothing
	end
	str *= ")"
	return Meta.parse(str)
end

@generated function possible_parent_cptids(id::CptID{N},psp::SVector{N,Int}) where {N}
#function possible_parent_cptids(::Type{Val{N}}) where {N}
	code = gen_possible_parent_cptids_code(Val{N})
	return quote
		SVector{$N,CptID{$N}}($code)
	end
end

#function possible_parent_cptids(id::CptID{N},psp::SVector{N,Int}) where {N}
#	ret = Vector{CptID{N}}(undef,N)
#	for i = 1:N
#		imulti = id[2]
#		ptidx = id[1]		
#		ptidx1d = parent_pt_idx(ptidx[i],imulti[i],psp[i])
#		imulti1d = imulti[i]-1
#		_imulti = _unroll_(imulti,imulti1d,Val{i})
#		_ptidx = _unroll_(ptidx,ptidx1d,Val{i})
#		ret[i] = (_ptidx,_imulti)
#	end
#	return SVector{N,CptID{N}}(ret)
#end


struct AncestorIterator1d{HCP<:AbstractHierarchicalCollocationPoint}
	cpt::HCP
	dim::Int
	stoplevel::Int
	function AncestorIterator1d(hcpt::HCP,dim::Int=1,stoplevel::Int=1) where {HCP<:AbstractHierarchicalCollocationPoint}
		return new{HCP}(hcpt,dim,stoplevel)
	end
end

cpt(iter::AncestorIterator1d) = iter.cpt
dim(iter::AncestorIterator1d) = iter.dim
stoplevel(iter::AncestorIterator1d) = iter.stoplevel
Base.eltype(iter::AncestorIterator1d{HCP}) where {HCP<:AbstractHierarchicalCollocationPoint} = HCP

function Base.length(iter::AncestorIterator1d{HCP}) where {HCP<:AbstractHierarchicalCollocationPoint}
	len = 0
	d = dim(iter)
	iter_cpt = cpt(iter)
	startlevel = i_multi(iter_cpt,d)
	stplvl = stoplevel(iter)
	par = parent(iter_cpt,d)
	while isvalid(par) && len < (startlevel-stplvl)
		len+=1
		par = parent(par,d)
	end
	return len
end



function Base.iterate(iter::AncestorIterator1d{HCP},state::HCP=cpt(iter)) where {HCP<:AbstractHierarchicalCollocationPoint}
	d = dim(iter)
	if i_multi(state,d) > stoplevel(iter)
		par = parent(state,d)
		if isvalid(par)
			return par,par
		else
			return nothing
		end
	else
		return nothing
	end
end


struct AncestorIterator{HCP<:AbstractHierarchicalCollocationPoint}
	cpt::HCP
	stoplevel::Int
	function AncestorIterator(cpt::HCP,stoplevel::Int=1) where {HCP<:AbstractHierarchicalCollocationPoint}
		return new{HCP}(cpt,stoplevel)
	end
end

cpt(iter::AncestorIterator) = iter.cpt
stoplevel(iter::AncestorIterator) = iter.stoplevel
Base.eltype(iter::AncestorIterator{HCP}) where {HCP<:AbstractHierarchicalCollocationPoint} = Set{HCP}

function Base.iterate(iter::AncestorIterator{HCP}) where {HCP<:AbstractHierarchicalCollocationPoint}
	i_cpt = cpt(iter)
	ancs = Set{HCP}()
	prnts = parents(i_cpt)
	for prntdim in prnts
		if isvalid(prntdim)
			push!(ancs, prntdim)
		end
	end
	return ancs,ancs
end

function Base.iterate(iter::AncestorIterator{HCP},state::Set{HCP}) where {HCP<:AbstractHierarchicalCollocationPoint}
	next_ancs = Set{HCP}()
	stplvl = stoplevel(iter)
	if !isempty(state) && stplvl < level(first(state))
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

function Base.length(iter::AncestorIterator{HCP}) where {HCP<:AbstractHierarchicalCollocationPoint}
	len = 0
	for elem in iter
		len+=1
	end
	return len
end


struct InterpolationIterator{N,CT,V<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	root_cpt::HCP
	x::V
	stoplevel::Int
	function InterpolationIterator(root::HCP, x::V, stoplevel::Int) where {N,CT,V<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
		return new{N,CT,V,CP,HCP}(root,x,stoplevel)
	end
end

Base.length(iter::InterpolationIterator) = begin; count = 0; for l in iter; count+=1; end; return count; end
stoplevel(iter::InterpolationIterator) = iter.stoplevel
get_root(iter::InterpolationIterator) = iter.root_cpt
coords(iter::InterpolationIterator) = iter.x
coord(iter::InterpolationIterator,dim::Int) = iter.x[dim]
Base.eltype(iter::InterpolationIterator{N,CT,V,CP,HCP}) where {N,CT,V<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}} = Set{HCP}

function Base.iterate(iter::InterpolationIterator{N,CT,V,CP,HCP}) where {N,CT,V<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	if stoplevel(iter) > 0
		root = get_root(iter)
		next_level_iter_cpts = Set{HCP}()
		push!(next_level_iter_cpts,root)
		return next_level_iter_cpts,next_level_iter_cpts
	else
		return nothing
	end
end

function Base.iterate(iter::InterpolationIterator{N,CT,V,CP,HCP},state::Set{HCP}) where {N,CT,V<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	next_level_iter_cpts = Set{HCP}()
	next_level_cpts = Set{HCP}()
	if stoplevel(iter) > level(first(state))
		for cpt in state
			if isrefined(cpt)
				for dim = 1:N
					push!(next_level_iter_cpts,next_interpolation_descendant(cpt,coord(iter,dim),dim))
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

#function Base.iterate(iter::InterpolationIterator{N,CT,V,CP,HCP}) where {N,CT,V<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#	if stoplevel(iter) > 0
#		root = get_root(iter)
#		next_level_iter_cpts = Set{HCP}()
#		state = Set{HCP}()
#		if isrefined(root) && stoplevel(iter) > 1
#			push!(state,root)
#		else
#			push!(next_level_iter_cpts,root)
#		end
#		if isempty(next_level_iter_cpts)
#			return iterate(iter, state)
#		else
#			return next_level_iter_cpts,state
#		end
#	else
#		return nothing
#	end
#end
#
#function Base.iterate(iter::InterpolationIterator{N,CT,V,CP,HCP},state::Set{HCP}) where {N,CT,V<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#	next_level_iter_cpts = Set{HCP}()
#	next_level_cpts = Set{HCP}()
#	if !isempty(state) && stoplevel(iter) > level(first(state))
#		for cpt in state
#			if isrefined(cpt)
#				for dim = 1:N
#					nli = next_interpolation_descendant(cpt,coord(iter,dim),dim)
#					if isrefined(nli) && stoplevel(iter) > level(nli)
#						push!(next_level_cpts, nli)
#					else
#						push!(next_level_iter_cpts,nli)
#					end
#				end
#			end
#		end
#	end
#	if isempty(next_level_cpts) && isempty(next_level_iter_cpts)
#		return nothing
#	else
#		if isempty(next_level_iter_cpts)
#			return iterate(iter, next_level_cpts)
#		else
#			return next_level_iter_cpts, next_level_cpts
#		end
#	end
#end

function next_interpolation_descendant(cpt::HCP,x::CT,dim::Int) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	chldrn = children(cpt,dim)
	vldch1 = isvalid(chldrn[1])
	vldch2 = isvalid(chldrn[2])
	crd = coord(cpt,dim)
	if vldch1 && vldch2
		if crd-x>zero(CT)
			#left child
			return cpt.children[dim][1]
		else
			#right child
			return cpt.children[dim][2]
		end
	elseif vldch1
		return chldrn[1]
	elseif vldch2
		return chldrn[2]
	else
		println("$cpt $x $dim")
		error()
	end
end

function Base.lastindex(iter::InterpolationIterator) 
	ret = iterate(iter)
	if isnothing(ret)
		error()
	end
	_lastset,__ = ret
	while !isnothing(ret)
	 	ret = iterate(iter,_lastset)
	 	if !isnothing(ret)
	 		_lastset,__ = ret
	 	end 
	 end
	 return _lastset
end