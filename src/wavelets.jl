multi_index(cpt) = cpt.i_multi

function wavelet_polyweight_exponents(cpt::CollocationPoint{N,CT}, numlevelsdown::Int, force_zero_mean::Bool=true) where {N,CT}
	lvl = level(cpt)
	@assert lvl > 1
	untillvl = lvl-numlevelsdown
	aci = AncestorIterator(cpt, untillvl)
	cpt_midx = multi_index(cpt)
	M = Vector{SVector{N,Int}}(undef, 0)
	Mm = Vector{SVector{N,Int}}(undef, 0)
	for aclvl in aci
		empty!(Mm)
		for ac in aclvl
			ac_midx = multi_index(ac)
			m = cpt_midx-ac_midx
			push!(Mm,m)
		end
		Mm = sort(Mm, rev=true)
		append!(M, Mm)
	end
	if force_zero_mean
		M[end] = zeros(SVector{N,Int})
	end
	return M
end

function wavelet_coeffs(cpt::CollocationPoint{N,CT}, maxp::Int, numlevelsdown::Int=1, force_zero_mean::Bool=true) where {N,CT}
	αs = wavelet_polyweight_exponents(cpt, numlevelsdown, force_zero_mean)
	untillvl = level(cpt)-numlevelsdown
	aci = AncestorIterator(cpt, untillvl)
	nvmoments = length(αs)
	A = zeros(nvmoments,nvmoments)
	b = zeros(nvmoments)
	for (imom,α) in enumerate(αs)
		ctr = 0
		for aclvl in aci
			for ac in aclvl
				ctr += 1
				A[imom, ctr] = integral_basis_fun(ac, maxp, α)
			end
		end
		b[imom] = integral_basis_fun(cpt, maxp, α)
	end
	#display(A)
	#display(b)
	return A\b
end

function wavelet(cpt::CollocationPoint{N,CT}, x::AV1, maxp::Int, numlevelsdown::Int=1, force_zero_mean::Bool=true, coeffs::AV2=wavelet_coeffs(cpt, maxp, numlevelsdown, force_zero_mean)) where {N,CT, AV1<:AbstractVector, AV2<:AbstractVector}
	lvl = level(cpt)
	untillvl = lvl-numlevelsdown
	aci = AncestorIterator(cpt, untillvl)
	res = basis_fun(cpt, x, maxp)
	ctr = 0
	for aclvl in aci
		for ac in aclvl
			ctr += 1
			res -= coeffs[ctr] * basis_fun(ac, x, maxp)
		end
	end
	return res
end

using SparseArrays, LinearAlgebra
function wavelet_to_hierarchical(asg::ASG{N,CT,Nv}, maxp::Int, numlevelsdown::Int=1) where {N,CT,Nv}
	(cpt_to_int, int_to_cpt) = dofmap(asg)
	npts = numpoints(asg)
	I = Vector{Int}(undef,0)
	J = Vector{Int}(undef,0)
	V = Vector{CT}(undef,0)
	push!(I,1)
	push!(J,1)
	push!(V,one(CT))
	wgts = zeros(numpoints(asg))
	#wgts[1] = weight(asg, first(values(asg.cpts[1])), x->x[1]^4, 3)
	for lvl in asg.cpts[2:end]
		for cpt in values(lvl)
			ind1 = cpt_to_int[getkey(cpt)]
			w = weight(asg, cpt, x->x[1]^4, 3)
			wgts[ind1] = w
			push!(I, ind1)
			push!(J, ind1)
			push!(V, one(CT))
			c = wavelet_coeffs(cpt, maxp, numlevelsdown)
			lvl = level(cpt)
			untillvl = lvl-numlevelsdown
			aci = AncestorIterator(cpt, untillvl)
			ctr = 0
			for lvl2 in aci
				for ac in lvl2
					ctr += 1
					ind2 = cpt_to_int[getkey(ac)]
					push!(I, ind2)
					push!(J, ind1)
					push!(V,-c[ctr])
				end
			end
		end
	end
	w2h = LinearAlgebra.UpperTriangular(SparseArrays.sparse(I,J,V))
	return cpt_to_int, int_to_cpt, w2h, wgts
end
