function numwaveletcoeffs(N::Int, numlevelsdown::Int)
	return (numlevelsdown+1)^N-1
end

function numlevelsdown_possible(hcpt::HCP, d::Int, numlevelsdown::Int=1) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	untillvl = level(hcpt,d)-numlevelsdown
	aci_d = AncestorIterator1d(hcpt, d, untillvl)
	return min(numlevelsdown, length(aci_d))
end

function correct_numlevelsdown_possible(hcpt::HCP, numlevelsdown_dims) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	if all(numlevelsdown_dims .== 0)
	#	return map(x->0,1:N)
		return numlevelsdown_dims
	end
	untillvl = maximum(map(d->level(hcpt,d)-numlevelsdown_dims[d],1:N))
	aci = AncestorIterator(hcpt, untillvl)
	#numlevelsdown_dim = ntuple(x->numlevelsdown,Val{N})
	len = sum(map(length,map(aclvl-> filter(ac->all(map((x,y)->x<=y,i_multi(hcpt)-i_multi(ac),numlevelsdown_dims)),aclvl), aci)))
	#map(aclvl-> filter(ac->all(map((x,y)->x<=y,i_multi(hcpt)-i_multi(ac),numlevelsdown_dims)),aclvl), aci)
	if reduce(*,map(x->x+1,numlevelsdown_dims))-1 == len
		return numlevelsdown_dims
	else
		println(hcpt)
		foreach(x->foreach(y->println(getkey(y)," ",i_multi(hcpt)-i_multi(y)," ",all(map((x,y)->x<=y,i_multi(hcpt)-i_multi(y),numlevelsdown_dims))),x),aci)
		foreach(x->println(x),hcpts_id_missing)
		error()
		_max = findmax(numlevelsdown_dims)[1]
		inds = findall(x->x==_max,numlevelsdown_dims)
		ind = inds[rand(1:length(inds))]
		@info "corrected numlevelsdown_dims $numlevelsdown_dims -> $(ntuple(n-> n == ind ? numlevelsdown_dims[n]-1 : numlevelsdown_dims[n] ,Val{N}))"
		new_numlevelsdown_dims = ntuple(n-> n == ind ? numlevelsdown_dims[n]-1 : numlevelsdown_dims[n] ,Val{N})
		
		return correct_numlevelsdown_possible(hcpt, new_numlevelsdown_dims)
	end
end


function get_missing_ten_prod_ids(asg, hcpt::HCP, numlevelsdown_dims) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	untillvl = maximum(map(d->level(hcpt,d)-numlevelsdown_dims[d],1:N))
	aci = AncestorIterator(hcpt, untillvl)
	hcpt_ids_available = Set{CptID}()
	hcpt_ids_neccessary = tens_prod_wave_cpts_ids(hcpt, 2)[2:end]
	filter!(x->all(map((x_,y_)->x_<=y_,i_multi(hcpt)-x[2],numlevelsdown_dims)),hcpt_ids_neccessary)
	for hcptset in aci
		for ac in hcptset
			if all(map((x,y)->x<=y,i_multi(hcpt)-i_multi(ac),numlevelsdown_dims))
				push!(hcpt_ids_available,getkey(ac))
			end
		end
	end
	hcpts_id_missing = setdiff(hcpt_ids_neccessary,hcpt_ids_available)
	return hcpts_id_missing
end

function collect_children(hcpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	allch = Set{HCP}()
	for chld in hcpt.children
		for ch in chld
			if isvalid(ch)
				push!(allch,ch)
				collect_children!(allch,ch)
			end
		end
	end
	return allch
end

function collect_children(allch::Set{HCP},hcpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	for chld in hcpt.children
		for ch in chld
			if isvalid(ch)
				push!(allch,ch)
				collect_children!(allch,ch)
			end
		end
	end
	return nothing
end

function force_numlevelsdown_possible_1d!(asg, hcpt::HCP, numlevelsdown, dim) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	for dim = 1:N
		untillvl = level(hcpt,dim)-numlevelsdown
		aci = AncestorIterator1d(hcpt, dim, untillvl)
		while length(aci) < numlevelsdown
	 		dimprnt_id = possible_parent_cptids(Val{dim}, getkey(hcpt), asg.pointSetProperties[dim])
	 		hsk = haskey(asg,dimprnt_id)
	 		while !hsk
	 			dimprnt_id = possible_parent_cptids(Val{dim}, dimprnt_id, asg.pointSetProperties[dim])
	 			hsk = haskey(asg,dimprnt_id)
	 		end
	 		npts = refine!(asg,asg[dimprnt_id])
		end
	end
end


function force_numlevelsdown_possible!(asg, hcpt::HCP, numlevelsdown_dims) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	
	untillvl = maximum(map(d->level(hcpt,d)-numlevelsdown_dims[d],1:N))
	hcpts_id_missing = get_missing_ten_prod_ids(asg, hcpt, numlevelsdown_dims)
	aci = AncestorIterator(hcpt, untillvl)
	new_points_created = false
	if !isempty(hcpts_id_missing)
		println("hcpt = $(getkey(hcpt))")
		println(hcpts_id_missing)
	end
	ctr = 0
	while !isempty(hcpts_id_missing) && ctr < 10
		ctr += 1
		bref = false
		for (lvl_down,hcptset) in enumerate(aci)
			npts_set = Set{CptID}()
			for ac in hcptset
				for idm in hcpts_id_missing
					if !isrefined(ac) && sum(idm[2]-i_multi(ac)) == lvl_down-sum(i_multi(hcpt)-idm[2]) && all(idm[2].>=i_multi(ac))
						println("refine $(getkey(ac))")
						npts = refine!(asg,ac)
						bref = true
						new_points_created = true
						for dimnpts in npts
							#for pt in dimnpts
								#setdiff!(hcpts_id_missing,map(x->getkey(x),dimnpts))
								union!(npts_set,map(x->getkey(x),dimnpts))
							#end
						end
					#elseif isrefined(ac) && sum(idm[2]-i_multi(ac)) == lvl_down-sum(i_multi(hcpt)-idm[2]) && all(idm[2].>=i_multi(ac))
					end
				end
			end
			if bref
				setdiff!(hcpts_id_missing,npts_set)
				break
			end
		end
	end

	aci = AncestorIterator(hcpt, untillvl)
	len = sum(map(length,map(aclvl-> filter(ac->all(map((x,y)->x<=y,i_multi(hcpt)-i_multi(ac),numlevelsdown_dims)),aclvl), aci)))
	if reduce(*,map(x->x+1,numlevelsdown_dims))-1 != len 
		println("numlevelsdown_dims = ",numlevelsdown_dims)
		println("len = ",len)
		println("hcpt = ",hcpt)
		println("ancs")
		foreach(x->foreach(y->println(getkey(y)," ",i_multi(hcpt)-i_multi(y)," ",all(map((x,y)->x<=y,i_multi(hcpt)-i_multi(y),numlevelsdown_dims))),x),aci)
		println("hcpt missing")
		println(hcpts_id_missing)
		error()
	end

	return new_points_created

end

function force_numlevelsdown_possible!(asg::AHSG{N,HCP},numlevelsdown::Int) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	b_something_refined = true
	ctr = 0
	while b_something_refined
		ctr+=1
		println("force_numlevelsdown_possible iteration $ctr")
		b_something_refined = false
		for hcpt in asg
			if all(i_multi(hcpt) .>= numlevelsdown+1)
				numlevelsdown_dims = ntuple(d->numlevelsdown,Val{N})
			else
				numlevelsdown_dims = ntuple(d->numlevelsdown_possible(hcpt, d, numlevelsdown),Val{N})
			end
			bsrtmp = force_numlevelsdown_possible!(asg, hcpt, numlevelsdown_dims)
			if bsrtmp 
				b_something_refined = true
			end
		end
	end
end

function tens_prod_wave_cpts_ids(cpt, numlevelsdown)
	cpt_idx, cpt_midx = getkey(cpt)
	N = length(cpt_idx)
	#nvmom = numlevelsdown_possible(cpt, numlevelsdown)
	#numlevelsdown = min(numlevelsdown,nvmom)
	tmp_idx = [Vector{Int}() for d = 1:N]
	tmp_midx = [Vector{Int}() for d = 1:N]



	for d = 1:N
		push!(tmp_idx[d], cpt_idx[d])
		push!(tmp_midx[d], cpt_midx[d])
		numlevelsdown = numlevelsdown_possible(cpt, d, numlevelsdown)
		untillvl = level(cpt,d)-numlevelsdown
		aci_d = AncestorIterator1d(cpt, d, untillvl)
		ac_ctr = 0
		for ac in aci_d
			ac_ctr += 1
			ac_idx, ac_midx = getkey(ac)
			push!(tmp_idx[d], ac_idx[d])
			push!(tmp_midx[d], ac_midx[d])
		end
	end
	idx = ndgrid(tmp_idx...)
	midx = ndgrid(tmp_midx...)
	res = Vector{CptID{N}}(length(midx[1]))
	for i = 1:length(midx[1])
		idxx = zeros(Int,N)
		midxx = zeros(Int,N)
		for d = 1:N
			idxx[d] = idx[d][i]
			midxx[d] = midx[d][i]
		end
		res[i] = (SVector{N}(idxx), SVector{N}(midxx))
	end
	return res
end

function tens_prod_wave_cpts(hcpt::HCP, numlevelsdown_dims) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	#res = Set{HCP}()
	res = Vector{HCP}(undef, reduce((x,y)->x*(y+1),numlevelsdown_dims,init=1)-1)
	untillvl_dim = [level(hcpt,d)-numlevelsdown_dims[d] for d in 1:N]
	untillvl = maximum(untillvl_dim)
	aci = AncestorIterator(hcpt, untillvl)
	ctr = 0
	for mac in aci
		for ac in mac
			if all(map((x,y)->x<=y,i_multi(hcpt)-i_multi(ac),numlevelsdown_dims))
				ctr += 1
				res[ctr] = ac
			end
		end
	end
	@assert ctr == length(res)
	inds = get_tprod_inds(hcpt, res, numlevelsdown_dims) 
	return res[inds]
end

function get_tprod_inds(hcpt::HCP, tpcpts::Vector{HCP}, numlevelsdown_dims) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	sze = map(x->x+1,numlevelsdown_dims)
	res = Vector{Int}(undef,length(tpcpts))
	ctr = 0
	for ac in tpcpts
		ctr += 1
		ii = (i_multi(hcpt)-i_multi(ac)+1)
		res[ctr] = sub2ind(sze,ii...)-1
	end
	return res
end

import StaticArrays: SMatrix
import LinearAlgebra
function wavelet_coeffs(asg, hcpt::HCP, maxp::Int, numlevelsdown::Int=1) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
				
		numlevelsdown_dims = ntuple(d->numlevelsdown_possible(hcpt, d, numlevelsdown),Val{N})
		#println(numlevelsdown_dims)
		

		#numlevelsdown_dims = correct_numlevelsdown_possible(hcpt, numlevelsdown_dims)
		
		#force_numlevelsdown_possible!(asg, hcpt, numlevelsdown_dims)
		
		#println(numlevelsdown_dims)
		#println()
		tpwcptids = tens_prod_wave_cpts(hcpt, numlevelsdown_dims)

		if all(numlevelsdown_dims .== 0)
			return zeros(0)
		end
	

		xd = [zeros(nvmom+1) for nvmom in numlevelsdown_dims]
		for d = 1:N
			xd[d][1] = 1.
			if numlevelsdown_dims[d] > 0
				xd[d][2:end] = -wavelet_coeffs_1d(asg, hcpt, d, maxp, numlevelsdown)
			end
		end
		dims = ntuple(i->numlevelsdown_dims[i]+1, Val{N}())
		nentrs = reduce(*,dims)
		res = zeros(nentrs)
		for i = 1:nentrs
			I = ind2sub(dims, i)
			res[i] = 1.
			for d = 1:N
				res[i] *= xd[d][I[d]]
			end
		end
		c = -res[2:end]
		
		if length(c) != length(tpwcptids)
			println(numlevelsdown_dims)
			foreach(x->println(getkey(x)),tpwcptids)
			error()
		end

		#display(c)
		#@info(integral_wavelet(asg, hcpt, maxp, numlevelsdown, zeros(SVector{N,Int}), c))
		return c
	
end



function wavelet_coeffs_1d(asg, hcpt::HCP, dim::Int, maxp::Int, numlevelsdown::Int=1) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	#nvmom = numlevelsdown_possible(hcpt, numlevelsdown)
	#numlevelsdown = min(numlevelsdown,nvmom)
	if numlevelsdown == 0
		return zeros(0)
	else
		numlevelsdown = numlevelsdown_possible(hcpt, dim, numlevelsdown)
		untillvl = level(hcpt,dim)-numlevelsdown
		aci_d = AncestorIterator1d(hcpt, dim, untillvl)
		A = Matrix{Float64}(undef, numlevelsdown, numlevelsdown)
		b = Vector{Float64}(undef, numlevelsdown)

		for imom = 1:numlevelsdown
			for (iac, ac) in enumerate(aci_d)
				A[imom, iac] = integral_basis_fun(ac, dim, maxp, true, imom-1)
			end
			b[imom] = integral_basis_fun(hcpt, dim, maxp, true, imom-1)
		end
		c = A\b
		for imom = 1:numlevelsdown
			res = integral_basis_fun(hcpt, dim, maxp, true, imom-1)
			for (iac, ac) in enumerate(aci_d)
				res -= c[iac] * integral_basis_fun(ac, dim, maxp, true, imom-1)
			end
			#warn("$(imom-1): $res")
			@assert isapprox(res, 0., atol=1e-12)
		end
		return c
	end
end






#
# function wavelet_coeffs!(asg, wcpt::WCP, maxp::Int, numlevelsdown::Int=1) where {N,CP,CT,RT,Nv,WCP<:WaveletCollocationPoint{N,CP,RT,CT,Nv}}
# 	w_coeffs = wavelet_coeffs(asg, wcpt, maxp, numlevelsdown)
# 	set_wavelet_coeffs!(wcpt,w_coeffs)
# end

function wavelet_coeffs!(wasg::SG, numlevelsdown::Int=1) where {N,CP,RT,CT,Nv,HCP<:WaveletCollocationPoint{N,CP,RT,CT,Nv},SG<:AdaptiveHierarchicalSparseGrid{N,HCP}}
	Maxp = maxporder(wasg)
	force_numlevelsdown_possible!(wasg,numlevelsdown)
	for wcpt in wasg
		if level(wcpt) > 1
			w_coeffs = wavelet_coeffs(wasg, wcpt, Maxp, numlevelsdown)
			set_wavelet_coeffs!(wcpt,w_coeffs)
		end

	end
end

function wavelet(asg, hcpt::HCP, x::AV1, maxp::Int, numlevelsdown::Int=1, coeffs::AV2=hcpt.wavelet_coeffs) where {N,CT,HCP<:WaveletCollocationPoint{N}, AV1<:AbstractVector{CT}, AV2<:AbstractVector}
	if all(map(x->isapprox(x,zero(CT),atol=1e-13),coeffs))
		return basis_fun(hcpt, x, maxp)
	else
		nvmom = numlevelsdown_possible(hcpt, numlevelsdown)
		numlevelsdown = min(numlevelsdown,nvmom)
		#res = 0.
		res = basis_fun(hcpt, x, maxp)
		cpts = map(x->asg[x], tens_prod_wave_cpts(hcpt, numlevelsdown))[2:end]
		for (icpt,cpt) in enumerate(cpts)
			res -= coeffs[icpt] * basis_fun(cpt, x, maxp)
		end
	end
	return res
end

function integral_wavelet(asg, hcpt::HCP, maxp::Int=3, numlevelsdown::Int=1,  weightorder::SVector{N,Int}=zeros(SVector{N,Int}), coeffs::AV=hcpt.wavelet_coeffs) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP},AV<:AbstractVector}
	if all(map(x->isapprox(x,zero(CT),atol=1e-13),coeffs))
		return integral_basis_fun(hcpt, maxp, weightorder)
	else
		res = integral_basis_fun(hcpt, maxp, weightorder)
		cpts = map(x->asg[x], tens_prod_wave_cpts(hcpt, numlevelsdown))[2:end]
		for (icpt,cpt) in enumerate(cpts)
			res -= coeffs[icpt] * integral_basis_fun(cpt, maxp, weightorder)
		end
	end
	return res
end

function dofmap(asg::SG) where {N, WCP<:WaveletCollocationPoint, SG<:AbstractHierarchicalSparseGrid{N,WCP}}
	ctr = 0
	npts = numpoints(asg)
	cpt_to_int = Dict{CptID{N}, Int}()
	sizehint!(cpt_to_int, npts)
	int_to_cpt = Vector{CptID{N}}(undef, npts)
	for lvl in asg.cpts
		for pdict in values(lvl)
			for cpt in values(pdict)
				cptk = getkey(cpt)
				ctr += 1
				cpt_to_int[cptk] = ctr
				int_to_cpt[ctr] = cptk
			end
		end
	end
	return (cpt_to_int, int_to_cpt)
end

function dofmap(asg::SG, lvl::Int) where {N, WCP<:WaveletCollocationPoint, SG<:AbstractHierarchicalSparseGrid{N,WCP}}
	ctr = 0
	npts = numpoints(asg)
	cpt_to_int = Dict{CptID{N}, Int}()
	sizehint!(cpt_to_int, npts)
	int_to_cpt = Vector{CptID{N}}(undef, npts)
	#for lvl in asg.cpts
		for pdict in values(asg.cpts[lvl])
			for cpt in values(pdict)
				cptk = getkey(cpt)
				ctr += 1
				cpt_to_int[cptk] = ctr
				int_to_cpt[ctr] = cptk
			end
		end
	#end
	return (cpt_to_int, int_to_cpt)
end

using SparseArrays, LinearAlgebra
function wavelet_to_hierarchical(asg::SG, maxp::Int, numlevelsdown::Int=1) where {N,CT,RT,CP<:AbstractCollocationPoint{N,CT},WCP<:WaveletCollocationPoint{N,CP,RT}, SG<:AbstractHierarchicalSparseGrid{N,WCP}}
	(cpt_to_int, int_to_cpt) = dofmap(asg)
	npts = numpoints(asg)
	I = Vector{Int}(undef,0)
	J = Vector{Int}(undef,0)
	V = Vector{CT}(undef,0)
	push!(I,1)
	push!(J,1)
	push!(V,one(CT))
	wgts = zeros(numpoints(asg))
	wgts[1] = scaling_weight(get_root(asg))
	iter_res = iterate(asg) # leave out root point
	if iter_res != nothing; iter_res = iterate(asg,iter_res[2]); end
	while iter_res != nothing
		cpt = iter_res[1]
		ind1 = cpt_to_int[getkey(cpt)]
		w = scaling_weight(cpt)
		wgts[ind1] = w
		push!(I, ind1)
		push!(J, ind1)
		push!(V, one(CT))
		c = cpt.wavelet_coeffs
		#println(c)
		#cpts = map(x->asg[x], tens_prod_wave_cpts(cpt, numlevelsdown))[2:end]
		#numlevelsdown_dims = ntuple(d->numlevelsdown_possible(cpt, d, numlevelsdown), Val{N}())
		numlevelsdown_dims = ntuple(d->numlevelsdown_possible(cpt, d, numlevelsdown),Val{N})
		numlevelsdown_dims = correct_numlevelsdown_possible(cpt, numlevelsdown_dims)
		cpts = tens_prod_wave_cpts(cpt, numlevelsdown_dims)
		#println(cpts)
		for (icpt,cpt2) in enumerate(cpts)
			ind2 = cpt_to_int[getkey(cpt2)]
			push!(I, ind2)
			push!(J, ind1)
			push!(V, -cpt.wavelet_coeffs[icpt])
		end
		iter_res = iterate(asg,iter_res[2])
	end
	w2h = LinearAlgebra.UpperTriangular(SparseArrays.sparse(I,J,V))
	#warn(det(w2h))
	#display(full(w2h))
	return cpt_to_int, int_to_cpt, w2h, wgts
end




function wavelet_weights!(asg::SG, f::F, numlevelsdown::Int=1) where {N,CP,RT,CT,Nv,HCP<:WaveletCollocationPoint{N,CP,RT,CT,Nv},SG<:AdaptiveHierarchicalSparseGrid{N,HCP},F<:Function}
	Maxp = maxporder(asg)
	nlvl = numlevels(asg)
	rp = get_root(asg)
	rp.wavelet_weight = rp.scaling_weight

		cpt_to_int, int_to_cpt, w2h, wgts = wavelet_to_hierarchical(asg, Maxp, numlevelsdown)

		wvlt_wghts = w2h\wgts

		for i = 1:length(wvlt_wghts)
			w = wvlt_wghts[i]
			k = int_to_cpt[i]
			cpt = asg[k]
			cpt.wavelet_weight = w
		end
	#end
	return nothing
end



#function wavelet(hcpt::HCP, x::AV1, maxp::Int, numlevelsdown::Int=1, force_zero_mean::Bool=true, coeffs::AV2=wavelet_coeffs(hcpt, maxp, numlevelsdown, force_zero_mean)) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}, AV1<:AbstractVector, AV2<:AbstractVector}
#	lvl = level(hcpt)
#	untillvl = lvl-numlevelsdown
#	aci = AncestorIterator(hcpt, untillvl)
#	res = basis_fun(hcpt, x, maxp)
#	ctr = 0
#	for aclvl in aci
#		for ac in aclvl
#			ctr += 1
#			res -= coeffs[ctr] * basis_fun(ac, x, maxp)
#		end
#	end
#	return res
#end
