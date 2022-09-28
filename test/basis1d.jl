function polyorder_v0(cpt::CollocationPoint{N,CT,Nv}, dim::Int, maxp::Int=level(cpt, dim)-1) where {N,CT<:Real,Nv}
	o = level(cpt, dim)-1
	return min(o, maxp)
end

function polybasis_v0(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::Float64, maxp::Int=3, level_1_constant::Bool=true) where {N, CT<:Real,Nv}
	# the polynomial basis as described in Pflüger's Diss at TUM
	lvl = level(cpt, dim)
	untillvl = lvl-maxp+2
	needsendpoints = lvl-1 > maxp
	if maxp > 1 || (lvl == 1 && level_1_constant)
		res = one(CT)
		x1 = coord(cpt, dim)
		_par = parent(cpt, dim)
		par = _par
		while isvalid(_par) && (!needsendpoints || level(_par, dim) >= untillvl)
			x0 = coord(_par, dim)
			res *= (x - x0) / (x1 - x0)
			par = _par
			_par = parent(_par, dim)
		end
		if needsendpoints
			if isvalid(par)
				interv = interval(cpt, dim)
			else
				interv = interval(par, dim)
			end
			x0 = interv[1]
			res *= (x - x0) / (x1 - x0)
			x0 = interv[2]
			res *= (x - x0) / (x1 - x0)
		end
		return res
	else
		return hat(cpt, dim, x)
	end
end

function polyorder_v1(cpt::CollocationPoint{N,CT,Nv}, dim::Int, maxp::Int, level_1_constant::Bool=true) where {N,CT,Nv}
	lvl = level(cpt,dim)
	if lvl == 1 && level_1_constant
		return 0
	elseif lvl == 1 && !level_1_constant
		return 1
	end
	if isrefined(cpt) # a non-leaf
		chldrn = children(cpt, dim)
		c1, c2 = chldrn[1], chldrn[2]
		# call recursive polyorder only on the refined children
		if isvalid(c1) && isvalid(c2)
			if isrefined(c1) && isrefined(c2)
				res = min(polyorder_v1(c1,dim,maxp), polyorder_v1(c2,dim,maxp))-1
			elseif isrefined(c1)
				res = polyorder_v1(c1,dim,maxp)-1
			elseif isrefined(c2)
				res = polyorder_v1(c2,dim,maxp)-1
			else
				res = min(polyorder_v1(c1,dim,maxp), polyorder_v1(c2,dim,maxp))-1
			end
		elseif isvalid(c1)
			if isrefined(c1)
				res = polyorder_v1(c1,dim,maxp)-1
			else
				res = polyorder_v1(parent(cpt,dim), dim, maxp)+1
			end
		else
			if isrefined(c2)
				res = polyorder_v1(c2,dim,maxp)-1
			else
				res = polyorder_v1(parent(cpt,dim), dim, maxp)+1
			end
		end
	else # a leaf
		par = parent(cpt,dim)
		parchldrn = children(par,dim)
		c1, c2 = parchldrn[1], parchldrn[2]
		# call recursive polyorder only on parents where only one child
		# is refined. in the case where cpt is a leaf which is not refined,
		# we won't end up in an infinite loop
		if isvalid(par)
			novalid_refined = (isvalid(c1) && !isrefined(c1) && isvalid(c2) && !isrefined(c2)) || (isvalid(c1) && !isrefined(c1) && !isvalid(c2)) || (isvalid(c2) && !isrefined(c2) && !isvalid(c1))
			if novalid_refined

				res = min(level(cpt,dim)-1, maxp)
				if res == 2
					println("found: $(level(cpt,dim)-1), $maxp")
				end
			else
				res = polyorder_v1(par,dim,maxp)+1
			end
		else # this must be the root point
			return 0
		end
	end
	return max(res, 1)
end

function polybasis_v1(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::Float64, maxp::Int=3, level_1_constant::Bool=true) where {N, CT<:Real,Nv}
	# the polynomial basis as described in Pflüger's Diss at TUM but pushed upwards
	lvl = level(cpt, dim)
	cptord = polyorder_v1(cpt, dim, maxp)
	if cptord != 1
		needsendpoints = lvl-1 > cptord
		nptsneeded = cptord
		#println(nptsneeded)
		nprtsneeded = nptsneeded
		if needsendpoints
			nprtsneeded -= 2
		end
		x1 = coord(cpt, dim)
		res = one(CT)
		if nprtsneeded > 0
			untillvl = lvl - nprtsneeded
			aci = AncestorIterator1d(cpt, dim, untillvl)
			@assert length(aci) == nprtsneeded "$(length(aci)) == $nprtsneeded"
			_ac, stat = iterate(aci)
			ac = _ac
			for i = 1:nprtsneeded
				ac = _ac
				x0 = coord(_ac, dim)
				res *= (x - x0) / (x1 - x0)
				if i < nprtsneeded
					_ac, stat = iterate(aci, stat)
				end
			end
		else
			ac = cpt
		end
		if needsendpoints
			interv = interval(ac, dim)
			x0 = interv[1]
			res *= (x - x0) / (x1 - x0)
			x0 = interv[2]
			res *= (x - x0) / (x1 - x0)
		end
	else
		res = hat(cpt, dim, x)
	end
	return res
end




function polybasis(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::Float64, maxp::Int=3, level_1_constant::Bool=true, variant::Int=1) where {N,CT,Nv}
	if variant == 0
		return polybasis_v0(cpt, dim, x, maxp, level_1_constant)
	elseif variant == 1
		return polybasis_v1(cpt, dim, x, maxp, level_1_constant)
	else
		error()
	end
end

function hat(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::CT) where {N, CT<:Real,Nv}
	interv = interval(cpt, dim)
	if interv[1] <= x <= interv[2]
		xcpt = coord(cpt, dim)
		if x <= xcpt && (xcpt - interv[1]) > eps()
			res = (x - interv[1]) / (xcpt - interv[1])
			return res
		else
			res = one(CT) - (x - xcpt) / (interv[2] - xcpt)
			return res
		end
	else
		return zero(CT)
	end
end




import FastGaussQuadrature: gausslegendre


function integral_basis_fun(cpt::CollocationPoint{N,CT,Nv}, dim::Int, maxp::Int=3, level_1_constant::Bool=true, weightorder::Int=0, variant::Int=1) where {N,CT<:Real,Nv}
	interv = cpt.interv[dim]
	weightorder = max(0, weightorder)
	if variant == 0
		cptord = polyorder_v0(cpt, dim, maxp)
	elseif variant == 1
		cptord = polyorder_v1(cpt, dim, maxp)
	else
		error()
	end
	o = cptord + weightorder
	nip = ceil(Int,(o+1)/2)
	xip_unit, wip = gausslegendre(nip)
	if cptord > 1
		Δx = interv[2]-interv[1]
		res = zero(CT)
		for (xu, w) in zip(xip_unit, wip)
			xip = (xu+one(xu)) / 2 * Δx + interv[1]
			bip = basis_fun(cpt, dim, xip, maxp, level_1_constant) * Δx
			res += xip^weightorder * bip*w
		end
	else
		Δx = cpt.coords[dim]-interv[1]
		res = zero(CT)
		for (xu, w) in zip(xip_unit, wip)
			xip = (xu+one(xu)) / 2 * Δx + interv[1]
			bip = basis_fun(cpt, dim, xip, maxp, level_1_constant) * Δx
			res += xip^weightorder * bip*w
			@assert !isnan(xip^weightorder)
			@assert !isnan(bip) "$xip"
		end
		Δx = interv[2]-cpt.coords[dim]
		for (xu, w) in zip(xip_unit, wip)
			xip = (xu+one(xu)) / 2 * Δx + cpt.coords[dim]
			bip = basis_fun(cpt, dim, xip, maxp, level_1_constant) * Δx
			res += xip^weightorder * bip*w
			@assert !isnan(xip^weightorder)
			@assert !isnan(bip)
		end

	end
	return res
end

function integral_basis_fun(cpt::CollocationPoint{N,CT,Nv}, maxp::Int=3, weightorder::SVector{N,Int}=zeros(SVector{N,Int})) where {N,CT<:Real,Nv}
	res = one(CT)
	for d = 1:N
		res *= integral_basis_fun(cpt, d, maxp, true, weightorder[d])
	end
	return res
end
