import LinearAlgebra: norm

function refine(::Type{Val{DIM}}, cpt::CP, psp::Int) where {DIM,N,CP<:AbstractCollocationPoint{N}}
	tmpch1,tmpch2 = invalid_collocation_point(CP),invalid_collocation_point(CP)
	((cl,l_interv,llvloff),(cr,r_interv,rlvloff)) = next_points(coord(cpt,DIM),interval(cpt,DIM),psp,i_multi(cpt,DIM))
	ptidx = isroot(cpt,DIM) ? 2 : next_level_pt_idx(pt_idx(cpt,DIM),psp)
	if isvalid_interv(l_interv)
		tmpch1 = CP(Val{DIM},cpt,cl,l_interv,ptidx-1)
	end
	if isvalid_interv(r_interv)
		tmpch2 = CP(Val{DIM},cpt,cr,r_interv,ptidx+1)
	end
	return SVector{2,CP}(tmpch1,tmpch2)
end

function refine!(::Type{Val{DIM}}, asg::AHSG{N,HCP}, hcpt::HCP) where {DIM,N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	psp = pointsetproperties(asg,DIM)
	tmpch1,tmpch2 = invalid_collocation_point(HCP),invalid_collocation_point(HCP)
	((cl,l_interv,llvloff),(cr,r_interv,rlvloff)) = next_points(coord(hcpt,DIM),interval(hcpt,DIM),psp,i_multi(hcpt,DIM))
	ptidx = isroot(hcpt,DIM) ? 2 : next_level_pt_idx(pt_idx(hcpt,DIM),psp)
	if isvalid_interv(l_interv)
		tmpch1 = HierarchicalCollocationPoint(Val{DIM},hcpt,cl,l_interv,ptidx-1)
		tmpch1 = push!(asg,tmpch1)
	end
	if isvalid_interv(r_interv)
		tmpch2 = HierarchicalCollocationPoint(Val{DIM},hcpt,cr,r_interv,ptidx+1)
		tmpch2 = push!(asg,tmpch2)
	end
	return SVector{2,HCP}(tmpch1,tmpch2)
end

function gen_refine_code(::Type{Val{N}}) where {N}
	str = "("
	for dim = 1:N
		str*="refine!(Val{$dim},asg,hcpt)"
		dim < N ? str*="," : nothing
	end
	str *= ")"
	return Meta.parse(str)
end

@generated function _unroll_refine_(asg::AHSG{N,HCP}, hcpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	refine_code = gen_refine_code(Val{N})
	return quote
		$refine_code
	end
end

function refine!(asg::AHSG{N,HCP}, hcpt::HCP) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	@assert !isrefined(hcpt) && isvalid(hcpt)
	nchildren = SVector{N,SVector{2,HCP}}(_unroll_refine_(asg,hcpt))
	setchildren!(hcpt,nchildren)
	return nchildren
end

"""
    generate_next_level!(asg::AHSG{N,HCP}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}

Generates all collocation point of the next hierarchical level.

# Constructor
- `asg::AHSG{N,HCP}`: sparse grid

"""	
function generate_next_level!(asg::AHSG{N,HCP}) where {N,HCP<:AbstractHierarchicalCollocationPoint{N}}
	nchildren = Set{HCP}()
	asglvl = asg.cpts[end]
	for dict_ptidx in values(asglvl)
		for cpt in values(dict_ptidx)
			if !isrefined(cpt)
				ncpts = refine!(asg,cpt)
				for cptdim in ncpts
					for hcpt in cptdim
						if isvalid(hcpt)
							push!(nchildren,hcpt)
						end
					end
				end
			end
		end
	end
	return nchildren
end

"""
    generate_next_level!(asg::AHSG{N,HCP}, tol::CT,maxlvl::Int) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}

Adaptively generate all collocation point of the next hierarchical level, where `norm(scaling_weight(cpt)) > tol && level(cpt)<maxlvl`.

# Constructor
- `asg::AHSG{N,HCP}`: sparse grid
- `tol::CT`: tolerance
- `maxlvl::Int`: maximum hierarchical level

"""	
function generate_next_level!(asg::AHSG{N,HCP}, tol::CT,maxlvl::Int) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	nchildren = Set{HCP}()
	asglvl = asg.cpts[end]
	for dict_ptidx in values(asglvl)
		for cpt in values(dict_ptidx)
			if !isrefined(cpt) && norm(scaling_weight(cpt)) > tol && level(cpt)<maxlvl
				ncpts = refine!(asg,cpt)
				for cptdim in ncpts
					for hcpt in cptdim
						if isvalid(hcpt)
							push!(nchildren,hcpt)
						end
					end
				end
			end
		end
	end
	return nchildren
end

function wavelet_refine!(asg::AHSG{N,WCP}, tol::CT, maxlvl::Int) where {N,CT,WCP<:WaveletCollocationPoint{}}
	torefine = Set{WCP}()
	nchildren = Set{WCP}()

	for (i,asglvl) in enumerate(asg.cpts)
		for dict_ptidx in values(asglvl)
			for cpt in values(dict_ptidx)
				if !isrefined(cpt) && abs(cpt.wavelet_weight) > tol && level(cpt) < maxlvl
					if abs(cpt.scaling_weight) <= tol
						#println("Level $i: $(getkey(cpt)) => Scal.w = $(abs(cpt.scaling_weight)) <= $tol, Wavl. W = $(abs(cpt.wavelet_weight)) > $tol")
					end
					push!(torefine,cpt)
				elseif !isrefined(cpt) && abs(cpt.scaling_weight) > tol
					#println("Level $i: $(getkey(cpt)) => Scal.w = $(abs(cpt.scaling_weight)) > $tol, Wavl. W = $(abs(cpt.wavelet_weight)) <= $tol")
				end
			end
		end
	end

	for cpt in torefine
		if !isrefined(cpt)
			ncpts = refine!(asg,cpt)
			for cptdim in ncpts
				for hcpt in cptdim
					if isvalid(hcpt)
						push!(nchildren,hcpt)
					end
				end
			end
		end
	end

	return nchildren
end
