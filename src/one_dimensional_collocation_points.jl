abstract type PointSetProperty end
abstract type ClosedPointSet <: PointSetProperty end
abstract type OpenPointSet <: PointSetProperty end
abstract type HalfOpenPointSet{LR} <: PointSetProperty end
isopen(psp::Type{HalfOpenPointSet{LR}}, cp::CT) where {LR ,CT<:Real} = LR == -1 ? (cp < zero(CT) ? true : false) : (LR == 1 ? (cp > zero(CT) ? true : false) : (error("HalfOpenPointSet{LR}: LR has to be either -1 (for left open) or 1 (for right open)")))
isopen(psp::Type{ClosedPointSet}, cp::CT) where {CT<:Real} = false
isopen(psp::Type{OpenPointSet}, cp::CT) where {CT<:Real} = true
isclosed(psp::Type{T}, cp::CT) where {CT<:Real,T<:PointSetProperty} = !isopen(psp,cp)

#1->CPS,2->OPS,3->LOPS,4->ROPS
isclosed(::Type{Val{1}}, cp::CT) where {CT<:Real} = isclosed(ClosedPointSet,cp)
isclosed(::Type{Val{2}}, cp::CT) where {CT<:Real} = isclosed(OpenPointSet,cp)
isclosed(::Type{Val{3}}, cp::CT) where {CT<:Real} = isclosed(HalfOpenPointSet{-1},cp)
isclosed(::Type{Val{4}}, cp::CT) where {CT<:Real} = isclosed(HalfOpenPointSet{1},cp)
isclosed(i::Int, cp::CT) where {CT<:Real} = i == 1 ? (isclosed(Val{1}, cp)) : ( i==2 ? (isclosed(Val{2}, cp)) : ( i==3 ? (isclosed(Val{3}, cp)) : (i==4 ? (isclosed(Val{4}, cp)) : (error()))) )

root_point(::Type{CT}) where {CT<:Real} = (zero(CT),(-one(CT),one(CT)))
left_point_level2(::Type{CT}) where {CT<:Real} = (-one(CT),(-one(CT),zero(CT)))
right_point_level2(::Type{CT}) where {CT<:Real} = (one(CT),(zero(CT),one(CT)))
next_left_point(pct::CT,p_interv::Union{Tuple{CT,CT},SVector{2,CT}}) where {CT<:Real} = ( p_interv[1] + (pct-p_interv[1])/2, (p_interv[1],pct) )
next_right_point(pct::CT,p_interv::Union{Tuple{CT,CT},SVector{2,CT}}) where {CT<:Real} = ( pct + (p_interv[2]-pct)/2, (pct,p_interv[2]) )
isvalid_interv(interv::Union{Tuple{CT,CT},SVector{2,CT}}) where {CT<:Real} = interv[1] < interv[2]

function next_points(pct::CT,p_interv::Union{Tuple{CT,CT},SVector{2,CT}}) where {CT<:Real}
	(lp,lp_interv) = next_left_point(pct,p_interv)
	(rp,rp_interv) = next_right_point(pct,p_interv)
	return (lp,lp_interv),(rp,rp_interv)
end

function next_points(pct::CT,p_interv::Union{Tuple{CT,CT},SVector{2,CT}},PSP::Int,lvl::Int) where {CT<:Real}
	llvloff = rlvloff = 0
	if lvl > 1 || (!isclosed(PSP,-one(CT)) && !isclosed(PSP,one(CT)))
		(lp,lp_interv),(rp,rp_interv) = next_points(pct,p_interv)
	elseif lvl == 1 && pct == zero(CT) #root point
		if isclosed(PSP,-one(CT)) #left
			(lp,lp_interv) = left_point_level2(CT)
			llvloff = 0
		else
			(lp,lp_interv) = next_left_point(pct,p_interv)
		end
		if isclosed(PSP,one(CT)) #right
			(rp,rp_interv) = right_point_level2(CT)
			rlvloff = 0
		else
			(rp,rp_interv) = next_right_point(pct,p_interv)
		end
	elseif lvl == 1
		if isclosed(PSP,pct)
			(lp,lp_interv) = left_point_level2(CT)
			(rp,rp_interv) = right_point_level2(CT)
		else
			(lp,lp_interv) = next_left_point(pct,p_interv)
			(rp,rp_interv) = next_right_point(pct,p_interv)
		end
	end
	return (lp,lp_interv,llvloff),(rp,rp_interv,rlvloff)
end

function next_level_pt_idx(pt_idx::Int, pt_lvl::Int, PSP::Int)
	#if pt_idx == 1
	#	return 2
	if pt_idx == pt_lvl == 1
		return 2
	elseif !isclosed(PSP,-1)
		pt_idx_nl = pt_idx*2
	else
		pt_idx_nl = pt_idx*2-1
	end
	return pt_idx_nl
end

function next_level_pt_idx(pt_idx::Int, pt_lvl::Int,  PSP::Int, leveloffset::Int)
	for i = 1:leveloffset
		pt_idx = next_level_pt_idx(pt_idx, pt_lvl, PSP)
	end
	return pt_idx
end
