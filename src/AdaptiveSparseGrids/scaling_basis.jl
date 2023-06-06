function polyorder(cpt::CP, dim::Int, maxp::Int=level(cpt, dim)-1) where {CP<:AbstractCollocationPoint}
#function polyorder_v0(cpt::CollocationPoint{N,CT,Nv}, dim::Int, maxp::Int=level(cpt, dim)-1) where {N,CT<:Real,Nv}
	o = level(cpt, dim)-1
	return min(o, maxp)
end
function polybasis_v0(cpt::HCP, dim::Int, x::CT, maxp::Int=3, level_1_constant::Bool=true) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#function polybasis_v0(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::Float64, maxp::Int=3, level_1_constant::Bool=true) where {N, CT<:Real,Nv}
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

function polyorder_v1(cpt::HCP, dim::Int, maxp::Int=3, level_1_constant::Bool=true) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#function polyorder_v1(cpt::CollocationPoint{N,CT,Nv}, dim::Int, maxp::Int, level_1_constant::Bool=true) where {N,CT,Nv}
	@assert isvalid(cpt)
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
				par = parent(cpt, dim)
				if isvalid(par)
					res = min(polyorder_v1(par, dim, maxp)+1, maxp)
				else
					res = 1
				end
			end
		else
			if isrefined(c2)
				res = polyorder_v1(c2,dim,maxp)-1
			else
				par = parent(cpt, dim)
				if isvalid(par)
					res = min(polyorder_v1(par, dim, maxp)+1, maxp)
				else
					res = 1
				end
			end
		end
	else # a leaf
		par = parent(cpt,dim)
		# call recursive polyorder only on parents where only one child
		# is refined. in the case where cpt is a leaf which is not refined,
		# we won't end up in an infinite loop
		if isvalid(par)
			parchldrn = children(par,dim)
			c1, c2 = parchldrn[1], parchldrn[2]
			novalid_refined = (isvalid(c1) && !isrefined(c1) && isvalid(c2) && !isrefined(c2)) || (isvalid(c1) && !isrefined(c1) && !isvalid(c2)) || (isvalid(c2) && !isrefined(c2) && !isvalid(c1))
			if novalid_refined

				res = min(level(cpt,dim)-1, maxp)
				if res == 2
					#println("found: $(level(cpt,dim)-1), $maxp")
				end
			else
				res = min(polyorder_v1(par,dim,maxp)+1, maxp)
			end
		else # doesn't have a parent in this direction due to inhomogeneity
			res = 1
			#error("$(level(cpt,dim)) $(isvalid(cpt)) $(getkey(cpt))  $(getkey(par))")
		end
	end
	return max(res, 1)
end

#function polybasis_v1(cpt::HCP, dim::Int, x::CT, maxp::Int=3, level_1_constant::Bool=true) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
##function polybasis_v1(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::Float64, maxp::Int=3, level_1_constant::Bool=true) where {N, CT<:Real,Nv}
#	# the polynomial basis as described in Pflüger's Diss at TUM but pushed upwards
#	lvl = level(cpt, dim)
#	idx = pt_idx(cpt,dim)
#	cptord = polyorder_v1(cpt, dim, maxp)
#	if cptord != 1
#		needsendpoints = lvl-1 > cptord
#		nptsneeded = cptord
#		#println(nptsneeded)
#		nprtsneeded = nptsneeded
#		if needsendpoints
#			nprtsneeded -= 2
#		end
#		x1 = coord(cpt, dim)
#		#println()
#		#println("x1: ", x1)
#		res = one(CT)
#		if nprtsneeded > 0
#			untillvl = lvl - nprtsneeded
#			#aci = AncestorIterator1d(cpt, dim, untillvl)
#			#@assert length(aci) == nprtsneeded "$(length(aci)) == $nprtsneeded"
#			@assert lvl > nprtsneeded "$(length(aci)) == $nprtsneeded"
#			#_ac, stat = iterate(aci)
#			_ac = parent_pt_idx_closed(idx, lvl)
#			ac = _ac
#			for i = 1:nprtsneeded
#				ac = _ac
#				x0 = CT(pt_coord_closed(_ac, lvl-i))
#				#x0 = coord(_ac, dim)
#				#println("x0: ", x0)
#				res *= (x - x0) / (x1 - x0)
#				if i < nprtsneeded
#					#_ac, stat = iterate(aci, stat)
#					_ac = parent_pt_idx_closed(_ac, lvl-i)
#				end
#			end
#		else
#			ac = idx
#		end
#		if needsendpoints
#			interv_lvl = lvl - nprtsneeded
#			interv = pt_interval_closed(ac, interv_lvl)
#			x0 = CT(interv[1])
#			if abs(x1-x0) < 10*eps()
#				println(interv_lvl," ", lvl, " ", cptord)
#			end
#			#if abs(x1-x0) > 10*eps()
#				#println("x0: ", x0)
#				res *= (x - x0) / (x1 - x0)
#			#else # set a double zero at the other
#			#	x0 = interv[2]
#			#	println("x0: ", x0)
#			#	res *= (x - x0) / (x1 - x0)
#			#end
#			x0 = CT(interv[2])
#			if abs(x1-x0) < 10*eps()
#				println(interv_lvl," ", lvl, " ", cptord)
#			end
#			#if abs(x1-x0) > 10*eps()
#				#println("x0: ", x0)
#				res *= (x - x0) / (x1 - x0)
#			#else # set a double zero at the other
#			#	x0 = interv[1]
#			#	println("x0: ", x0)
#			#	res *= (x - x0) / (x1 - x0)
#			#end
#		end
#	else
#		res = hat(cpt, dim, x)
#	end
#	@assert !isinf(res) && !isnan(res)
#	return res
#end





function polybasis(cpt::HCP, dim::Int, x::CT, maxp::Int=3, level_1_constant::Bool=true, variant::Int=1) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#function polybasis(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::Float64, maxp::Int=3, level_1_constant::Bool=true, variant::Int=1) where {N,CT,Nv}
	if variant == 0
		return polybasis_v0(cpt, dim, x, maxp, level_1_constant)
	elseif variant == 1
		return polybasis_v1(cpt, dim, x, maxp, level_1_constant)
	else
		error()
	end
end

function hat(cpt::HCP, dim::Int, x::CT) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#function hat(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::CT) where {N, CT<:Real,Nv}
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


function derivative_hat(cpt::HCP, dim::Int, x::CT) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#function hat(cpt::CollocationPoint{N,CT,Nv}, dim::Int, x::CT) where {N, CT<:Real,Nv}
	interv = interval(cpt, dim)
	if interv[1] <= x <= interv[2]
		xcpt = coord(cpt, dim)
		if x <= xcpt && (xcpt - interv[1]) > eps()
			res = (1 - interv[1]) / (xcpt - interv[1])
			return res
		else
			res = - 1.0 / (interv[2] - xcpt)
			return res
		end
	else
		return zero(CT)
	end
end



import FastGaussQuadrature: gausslegendre

#function integral_basis_fun(_cpt::HCP, _dim::Int, interv::SVector{2,F}, xip_unit, wip, maxp::Int=3, level_1_constant::Bool=true, weightorder::Int=0) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP},F}
#	res = zero(CT)
#	if isvalid_interv(interv)
#		Δx = interv[2]-interv[1]
#		for (xu, w) in zip(xip_unit, wip)
#			xip = (xu+one(xu)) / 2 * Δx + interv[1]
#			bip = basis_fun(_cpt, _dim, xip, maxp, level_1_constant) #* Δx
#			res += xip^weightorder * bip * w * Δx/2.0
#			@assert !isnan(xip^weightorder) && !isinf(xip^weightorder)
#			@assert !isnan(bip) && !isinf(bip)
#		end
#	end
#	return res
#end

function inner_basis_fun(icpt::HCP, jcpt::HCP, _dim::Int, interv::SVector{2,F}, xip_unit, wip, maxp::Int=3, level_1_constant::Bool=true) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP},F}

	res = zero(CT)
	if isvalid_interv(interv)
		Δx = interv[2]-interv[1]
		for (xu, w) in zip(xip_unit, wip)
			xip = (xu+one(xu)) / 2 * Δx + interv[1]
			ibip = basis_fun(icpt, _dim, xip, maxp, level_1_constant) #* Δx
			jbip = basis_fun(jcpt, _dim, xip, maxp, level_1_constant) #* Δx
			res += ibip * jbip * w * Δx/2.0
			@assert !isnan(ibip) && !isinf(ibip)
			@assert !isnan(jbip) && !isinf(jbip)
		end
	end
	return res
end

function split_interval(a::SVector{2,F}) where {F}
	c = (a[2]-a[1])/2 + a[1]
	oa = SVector{2,F}((a[1],c))
	ob = SVector{2,F}((c,a[2]))
	return (oa,ob)
end

#function integral_basis_fun(_cpt::HCP, _dim::Int, maxp::Int=3, level_1_constant::Bool=true, weightorder::Int=0, variant::Int=1) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
##function integral_basis_fun(cpt::CollocationPoint{N,CT,Nv}, _dim::Int, maxp::Int=3, level_1_constant::Bool=true, weightorder::Int=0, variant::Int=1) where {N,CT<:Real,Nv}
#	#interv = cpt.interv[_dim]
#	interv = interval(_cpt,_dim)
#	weightorder = max(0, weightorder)
#	if variant == 0
#		cptord = polyorder_v0(_cpt, _dim, maxp)
#	elseif variant == 1
#		cptord = polyorder_v1(_cpt, _dim, maxp)
#	else
#		error()
#	end
#	o = cptord + weightorder
#	nip = ceil(Int,(o+1)/2)
#	xip_unit, wip = gausslegendre(nip)
#	if cptord > 1
#		res = integral_basis_fun(_cpt, _dim, interv, xip_unit, wip, maxp, level_1_constant, weightorder)
#	else
#		(intva, intvb) = split_interval(interv)
#		res = integral_basis_fun(_cpt, _dim, intva, xip_unit, wip, maxp, level_1_constant, weightorder)
#		res += integral_basis_fun(_cpt, _dim, intvb, xip_unit, wip, maxp, level_1_constant, weightorder)
#	end
#	@assert !isnan(res) && !isinf(res) "$res"
#	return res
#end

function intersect_intervals(a::SVector{2,F}, b::SVector{2,F}) where {F}
	aisv = isvalid_interv(a)
	bisv = isvalid_interv(b)
	isdisjoint = a[1] >= b[2] || b[1] >= a[2]
	if aisv && bisv && !isdisjoint
		if a[1] <= b[1] && b[2] <= a[2] # b is contained in a
			ob = SVector{2,F}((b[1],b[2]))
		elseif b[1] <= a[1] && a[2] <= b[2] # a is contained in b
			ob = SVector{2,F}((a[1],a[2]))
		elseif b[1] < a[1] # b is left of a
			ob = SVector{2,F}((a[1],b[2]))
		elseif a[1] < b[1] # a is left of b
			ob = SVector{2,F}((b[1],a[2]))
		else
			error("$a, $b")
		end
	else
		ob = SVector{2,F}((1.,0.))
	end
	return ob
end


#function inner_basis_fun(icpt::HCP, jcpt::HCP, _dim::Int, maxp::Int=3, level_1_constant::Bool=true, variant::Int=1) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
##function integral_basis_fun(cpt::CollocationPoint{N,CT,Nv}, _dim::Int, maxp::Int=3, level_1_constant::Bool=true, weightorder::Int=0, variant::Int=1) where {N,CT<:Real,Nv}
#	#interv = cpt.interv[_dim]
#	i_interv = interval(icpt,_dim)
#	j_interv = interval(jcpt,_dim)
#	if variant == 0
#		icptord = polyorder_v0(icpt, _dim, maxp)
#		jcptord = polyorder_v0(jcpt, _dim, maxp)
#	elseif variant == 1
#		icptord = polyorder_v1(icpt, _dim, maxp)
#		jcptord = polyorder_v1(jcpt, _dim, maxp)
#	else
#		error()
#	end
#	o = icptord + jcptord
#	nip = ceil(Int,(o+1)/2)
#	xip_unit, wip = gausslegendre(nip)
#	res = zero(CT)
#	if icptord != 1 && jcptord != 1
#		int_intv = intersect_intervals(i_interv, j_interv)
#		res += inner_basis_fun(icpt, jcpt, _dim, int_intv, xip_unit, wip, maxp, level_1_constant)
#	elseif icptord == 1 && jcptord != 1
#		(ia,ib) = split_interval(i_interv)
#		aint_intv = intersect_intervals(ia, j_interv)
#		bint_intv = intersect_intervals(ib, j_interv)
#		int_intv_i = (intersect_intervals(ia, j_interv), intersect_intervals(ib, j_interv))
#		for intv in int_intv_i
#			res += inner_basis_fun(icpt, jcpt, _dim, intv, xip_unit, wip, maxp, level_1_constant)
#		end
#	elseif icptord != 1 && jcptord == 1
#		(ja,jb) = split_interval(j_interv)
#		aint_intv = intersect_intervals(i_interv, ja)
#		bint_intv = intersect_intervals(i_interv, jb)
#		int_intv_j = (intersect_intervals(i_interv, ja), intersect_intervals(i_interv, jb))
#		for intv in int_intv_j
#			res += inner_basis_fun(icpt, jcpt, _dim, intv, xip_unit, wip, maxp, level_1_constant)
#		end
#	elseif icptord == 1 && jcptord == 1
#		(ia,ib) = split_interval(i_interv)
#		(ja,jb) = split_interval(j_interv)
#		int_intv_ij = 	(intersect_intervals(ia, ja), intersect_intervals(ib, ja),
#						 intersect_intervals(ia, jb), intersect_intervals(ib, jb))
#		for intv in int_intv_ij
#			res += inner_basis_fun(icpt, jcpt, _dim, intv, xip_unit, wip, maxp, level_1_constant)
#		end
#	else
#		error()
#	end
#	return res
#end

#function basis_fun(hcpt::HCP, dim::Int, x::CT, maxp::Int=3, level_1_constant::Bool=true) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#	interv = interval(hcpt,dim)
#	if interv[1] <= x <= interv[2]
#		#if cpt.level > maxp
#		#	return zero(CT)
#		#else
#		#@assert norm(lagrange_basis_polynomial(cpt, dim, x, maxp, level_1_constant) - polybasis(cpt, dim, x, maxp, level_1_constant)) <= 100*eps()
#			return polybasis(hcpt, dim, x, maxp, level_1_constant)
#		#end
#	else
#		return zero(CT)
#	end
#end

function basis_fun(hcpt::HCP, x::AV, maxp::Int=3, level_1_constant::Bool=true) where {N,CT,AV<:AbstractVector{CT},CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	@assert length(x) == N
	@inbounds res = basis_fun(hcpt, 1, x[1], maxp, level_1_constant)
	@inbounds for dim = 2:N
		res *= basis_fun(hcpt, dim, x[dim], maxp, level_1_constant)
	end
	return res
end

#function integral_basis_fun(hcpt::HCP, maxp::Int=3, weightorder::SVector{N,Int}=zeros(SVector{N,Int})) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#	res = one(CT)
#	for d = 1:N
#		res *= integral_basis_fun(hcpt, d, maxp, true, weightorder[d])
#	end
#	return res
#end

#function inner_basis_fun(icpt::HCP, jcpt::HCP, maxp::Int=3) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
#	res = one(CT)
#	for d = 1:N
#		res *= inner_basis_fun(icpt, jcpt, d, maxp, true)
#	end
#	return res
#end

function basis_fun(hcpt::HCP, dim::Int, x::CT, maxp::Int=3, level_1_constant::Bool=true) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	interv = interval(hcpt,dim)
	if interv[1] <= x <= interv[2]
		#if cpt.level > maxp
		#	return zero(CT)
		#else
		#@assert norm(lagrange_basis_polynomial(cpt, dim, x, maxp, level_1_constant) - polybasis(cpt, dim, x, maxp, level_1_constant)) <= 100*eps()
			#return polybasis(hcpt, dim, x, maxp, level_1_constant)
		#end
		if i_multi(hcpt)[dim] == 1
			return one(CT)
		else
			return hat(hcpt, dim, x)
		end
	else
		return zero(CT)
	end
end

function integral_basis_fun(hcpt::HCP) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	return integral_basis_fun(cpt(hcpt))
end

function integral_basis_fun(hcpt::HCP,_dim::Int) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	return integral_basis_fun(cpt(hcpt),_dim)
end

function integral_basis_fun(cpt::CollocationPoint{N,CT}) where {N,CT}
	a = one(CT)
	for i = 1:N
		a *= integral_basis_fun(cpt,i)
	end
	return a
end

function integral_basis_fun(cpt::CollocationPoint{N,CT}, _dim::Int) where {N,CT}
	if i_multi(cpt)[_dim] == 1
		return CT(2.0)
	else
		i_interv = interval(cpt,_dim)
		return abs(i_interv[1]-i_interv[2])/CT(2.0)
	end
end


function derivative_basis_fun(hcpt::HCP,_dim::Int,x::CT) where {N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	return derivative_hat(hcpt, _dim, x)
end