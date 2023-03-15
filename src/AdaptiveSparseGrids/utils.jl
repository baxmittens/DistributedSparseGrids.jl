"""
	const AHSG{N,HCP} = AdaptiveHierarchicalSparseGrid{N,HCP}

Shortcut for AdaptiveHierarchicalSparseGrid

"""	
const AHSG{N,HCP} = AdaptiveHierarchicalSparseGrid{N,HCP}

get_root(asg::AHSG) = first(values(first(values(asg.cpts[1]))))

numlevels(cpts::Vector{PointDict{N,HCP}}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}} = length(cpts)
numlevels(asg::AHSG) = numlevels(asg.cpts)

pointsetproperties(asg::AHSG) = asg.pointSetProperties
pointsetproperties(asg::AHSG,dim::Int) = asg.pointSetProperties[dim]

maxporder(asg::AHSG) = asg.maxp

function Base.getindex(hcpts::PointDict{N, HCP}, ind::CptID{N}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	return hcpts[ind[2]][ind[1]]
end

function Base.getindex(hcpts::Vector{PointDict{N, HCP}}, ind::CptID{N}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	lvl = level(ind)
	return getindex(hcpts[lvl], ind)
end

function Base.getindex(asg::AHSG{N,CP}, ind::CptID{N}) where {N,CP<:AbstractHierarchicalCollocationPoint{N}}
	return getindex(asg.cpts, ind)
end

function Base.push!(ptidx_dict::Dict{SVector{N,Int},HCP},hcpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	ptidx = pt_idx(hcpt)
	if haskey(ptidx_dict, ptidx)
		retcpt = ptidx_dict[ptidx]
		pind = -1
		for dim = 1:N; if isvalid(parent(hcpt,dim)); pind = dim; break; end; end
		setparents!(retcpt, _unroll_dyn_dim_(parents(retcpt),parent(hcpt,pind),pind))
		return retcpt
	else
		ptidx_dict[ptidx] = hcpt
		return hcpt
	end
end

function Base.push!(ptdict::PointDict{N, HCP},hcpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	imul = i_multi(hcpt)
	if haskey(ptdict, imul)
		ptidx_dict = ptdict[imul]
		return push!(ptidx_dict, hcpt)
	else
		ptidx_dict = Dict{SVector{N,Int},HCP}()
		ptdict[imul] = ptidx_dict
		return push!(ptidx_dict, hcpt)
	end
end

function Base.push!(hcpts::Vector{PointDict{N, HCP}},hcpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	lvl = level(hcpt)
	nlvl = numlevels(hcpts)
	if lvl > nlvl
		@assert lvl == nlvl+1
		push!(hcpts,PointDict{N, HCP}())
	end
	return push!(hcpts[lvl],hcpt)
end

function Base.push!(asg::AHSG{N,HCP},_hcpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	hcpt = push!(asg.cpts,_hcpt)
	if level(hcpt) > 2
		ppids = possible_parent_cptids(getkey(hcpt),asg.pointSetProperties)
		
		for (dim,pid) in enumerate(ppids)
			if !isvalid(parent(hcpt,dim)) && haskey(asg,pid)
				prnt = asg[pid]
				setparents!(hcpt, _unroll_dyn_dim_(parents(hcpt), prnt, dim))
				lr = pt_idx(hcpt,dim) < next_level_pt_idx(pt_idx(prnt,dim),i_multi(prnt,dim),asg.pointSetProperties[dim]) ? 1 : 2
				setchild!(prnt, Val{dim}, Val{lr} ,hcpt)
				#println("dim=$dim: parent $pid set to $(getkey(hcpt))")
			end
		end
		nchildern = Vector{SVector{2,HCP}}()
		chg = false
		for i = 1:N
			possible_childs = refine(Val{i}, cpt(hcpt), asg.pointSetProperties[i])
			tmpch1,tmpch2 = children(hcpt,i)
			if getkey(possible_childs[1]) != getkey(tmpch1) && isvalid(possible_childs[1]) && haskey(asg,getkey(possible_childs[1]))
				tmpch1 = asg[getkey(possible_childs[1])]
				setparents!(tmpch1, _unroll_dyn_dim_(parents(tmpch1), hcpt, i))
				chg = true
			end
			if getkey(possible_childs[2]) != getkey(tmpch2) && isvalid(possible_childs[2]) && haskey(asg,getkey(possible_childs[2]))
				tmpch2 = asg[getkey(possible_childs[2])]
				setparents!(tmpch2, _unroll_dyn_dim_(parents(tmpch2), hcpt, i))
				chg = true
			end
			push!(nchildern,SVector{2,HCP}(tmpch1,tmpch2))
		end
		if chg
			#println("Update children")
			setchildren!(hcpt, SVector{N,SVector{2,HCP}}(nchildern))
		end
	end
	return hcpt
end

function Base.haskey(ptdict::PointDict{N,HCP},id::CptID{N}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	imul = id[2]
	if haskey(ptdict,imul)
		return haskey(ptdict[imul],id[1])
	else
		return false
	end
end

function Base.haskey(asg::AHSG{N,HCP},id::CptID{N}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	lvl = level(id)
	if any(map(x->x==0,id[2])) || lvl > numlevels(asg)
		return false
	else
		return haskey(asg.cpts[lvl],id)
	end
end

function Base.iterate(hcpts::PointDict{N,HCP}) where {N,HCP<: AbstractHierarchicalCollocationPoint{N}}

end

Base.eltype(asg::AHSG{N,HCP}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}} = HCP
Base.length(asg::AHSG) = numpoints(asg)
function Base.iterate(asg::AHSG{N,HCP}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	res_level_dict,state_level = iterate(asg.cpts)
	res_level = values(res_level_dict)
	res_imulti_dict,state_imulti = iterate(res_level)
	res_imulti = values(res_imulti_dict)
	res_ptidx,state_ptidx = iterate(res_imulti)
	return res_ptidx,(res_level,res_imulti,state_level,state_imulti,state_ptidx)
end
function Base.iterate(asg::AHSG{N,HCP},state) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	res_level,res_imulti,state_level,state_imulti,state_ptidx = state
	iter_ptidx_res = iterate(res_imulti,state_ptidx)
	if iter_ptidx_res != nothing
		res_ptidx,state_ptidx = iter_ptidx_res
		return res_ptidx,(res_level,res_imulti,state_level,state_imulti,state_ptidx)
	else
		iter_imulti_res = iterate(res_level,state_imulti)
		if iter_imulti_res != nothing
			res_imulti_dict,state_imulti = iter_imulti_res
			res_imulti = values(res_imulti_dict)
			res_ptidx,state_ptidx = iterate(res_imulti)
			return res_ptidx,(res_level,res_imulti,state_level,state_imulti,state_ptidx)
		else
			iter_level_res = iterate(asg.cpts,state_level)
			if iter_level_res != nothing
				res_level_dict,state_level = iter_level_res
				res_level = values(res_level_dict)
				res_imulti_dict,state_imulti = iterate(res_level)
				res_imulti = values(res_imulti_dict)
				res_ptidx,state_ptidx = iterate(res_imulti)
				return res_ptidx,(res_level,res_imulti,state_level,state_imulti,state_ptidx)
			end
		end
	end
	return nothing
end

function numpoints(asg::SG) where {SG<:AbstractHierarchicalSparseGrid}
	res = 0
	for level in asg.cpts
		for imultidict in values(level)
			res += length(imultidict)
		end
	end
	return res
end

function numpoints(asg::SG,lvl::Int) where {SG<:AbstractHierarchicalSparseGrid}
	res = 0
	#for level in asg.cpts
		for imultidict in values(asg.cpts[lvl])
			res += length(imultidict)
		end
	#end
	return res
end

function InterpolationIterator(asg::AHSG{N,HCP}, x::V, stoplevel::Int=numlevels(asg)) where {N,CT,V<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	return InterpolationIterator(get_root(asg),x,stoplevel)
end

function average_scaling_weight(asg::AHSG{N,HCP},lvl::Int) where {N,CT,V<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
    avg_val = similar(scaling_weight(first(asg)))
    tmp = similar(scaling_weight(first(asg)))
    n = 0
    for cpts in asg
        if level(cpts) == lvl
            n += 1
            fill!(tmp, 0.0)
            add!(tmp, scaling_weight(cpts))
            pow!(tmp, 2.0)
            pow!(tmp, 0.5)
            add!(avg_val, tmp)
        end
    end
    #pow!(avg_val, 0.5)
    mul!(avg_val, 1.0/n)
    return avg_val
end