
include("./CollocationPoints/one_dimensional_collocation_points.jl")
include("./CollocationPoints/tuple_generation.jl")

abstract type AbstractCollocationPoint{N,CT<:Real} end
abstract type AbstractHierarchicalCollocationPoint{N,CP<:AbstractCollocationPoint,RT} end


"""
	CollocationPoint{N,CT} <: AbstractCollocationPoint{N,CT}	
	
A collocation point.

# Fields

`i_multi::SVector{N,Int}` : The i-th item of this Vector represents to the level of the collocation point in the i-th dimension.
`pt_idx::SVector{N,Int}` : The i-th item of this Vector represents to the point index of the collocation point in the i-th dimension.
`coords::SVector{N,CT}` : Coordinates of the collocation point.
`interv::SVector{N,SVector{2,CT}}` : The i-th item of the Vector defines the non-zero interval of the associated basis function in the i-th dimension.
"""
struct CollocationPoint{N,CT} <: AbstractCollocationPoint{N,CT}	
	i_multi::SVector{N,Int}
	pt_idx::SVector{N,Int}
	coords::SVector{N,CT}
	interv::SVector{N,SVector{2,CT}}
end

include("./CollocationPoints/utils_CollocationPoint.jl")

function CollocationPoint{N,CT}(::Type{Val{DIM}}, parent::CollocationPoint{N,CT}, coord::CT, cp_interv::Tuple{CT,CT}, ptidx::Int) where {DIM,N,CT<:Real}
	return CollocationPoint{N,CT}(
		_unroll_(parent.i_multi,parent.i_multi[DIM]+1,Val{DIM}),
		_unroll_(parent.pt_idx,ptidx,Val{DIM}),
		_unroll_(parent.coords,coord,Val{DIM}),
		_unroll_(parent.interv,SVector{2,CT}(cp_interv),Val{DIM})
	)
end

@generated function invalid_collocation_point(::Type{CollocationPoint{N,CT}}) where {N,CT<:Real}
	sero = zero(CT)
	cpt = CollocationPoint{N,CT}(
		_unroll_(SVector{N,Int},-1),
		_unroll_(SVector{N,Int},-1),
		_unroll_(SVector{N,CT},sero),
		_unroll_(SVector{N,SVector{2,CT}},SVector{2,CT}(sero,sero))
		)
	return quote
		$cpt
	end
end

"""
    root_point(::Type{CollocationPoint{N,CT}}) where {N,CT<:Real}

generates a root point.

# Constructor
- `::Type{CollocationPoint{N,CT}}`: a colloction point of dimension `N` and collocation type `CT`.

"""	
@generated function root_point(::Type{CollocationPoint{N,CT}}) where {N,CT<:Real}
	rp,rp_interv = root_point(CT)
	cpt = CollocationPoint{N,CT}(
		_unroll_(SVector{N,Int},1),
		_unroll_(SVector{N,Int},1),
		_unroll_(SVector{N,CT},rp),
		_unroll_(SVector{N,SVector{2,CT}},SVector{2,CT}(rp_interv)),
		)
	return quote
		$cpt
	end
end


#if  ~isdefined(:Main,:HierarchicalCollocationPoint) #top level hack

"""
	HierarchicalCollocationPoint{N,CP,RT} <: AbstractHierarchicalCollocationPoint{N,CP,RT}	
	
	A collocation point.

	# Fields

	`cpt::CP` : [`DistributedSparseGrids.CollocationPoint`](@ref)
	`children::SVector{N,SVector{2,HierarchicalCollocationPoint{N,CP,RT}}}` : Container for at most two children per dimension
	`parents::SVector{N,HierarchicalCollocationPoint{N,CP,RT}}` : Container for parents 
	`fval::RT` : Function Value
	`scaling_weight::RT` : Scaling Weight (function value minus l-1 level interpolator both at cpt.coords)
"""
mutable struct HierarchicalCollocationPoint{N,CP,RT} <: AbstractHierarchicalCollocationPoint{N,CP,RT}
	cpt::CP
	children::SVector{N,SVector{2,HierarchicalCollocationPoint{N,CP,RT}}}
	parents::SVector{N,HierarchicalCollocationPoint{N,CP,RT}}
	fval::RT
	scaling_weight::RT
	function HierarchicalCollocationPoint{N,CP,RT}(cpt::CP) where {N,CP<:AbstractCollocationPoint{N},RT}
		return new{N,CP,RT}(cpt)
	end
	function HierarchicalCollocationPoint{N,CP,RT}(cpt::CP,children::SVector{N,SVector{2,HierarchicalCollocationPoint{N,CP,RT}}},parents::SVector{N,HierarchicalCollocationPoint{N,CP,RT}}) where {N,CP<:AbstractCollocationPoint{N},RT}
		return new{N,CP,RT}(cpt,children,parents)
	end
end
#end

if  ~isdefined(:Main,:WaveletCollocationPoint) #top level hack
mutable struct WaveletCollocationPoint{N,CP,RT,CT,Nv} <: AbstractHierarchicalCollocationPoint{N,CP,RT} 
	cpt::CP
	children::SVector{N,SVector{2,WaveletCollocationPoint{N,CP,RT,CT,Nv}}}
	parents::SVector{N,WaveletCollocationPoint{N,CP,RT,CT,Nv}}
	scaling_weight::RT
	wavelet_coeffs::SVector{Nv,CT}
	wavelet_weight::RT
	function WaveletCollocationPoint{N,CP,RT,CT,Nv}(cpt::CP) where {N,CT,CP<:AbstractCollocationPoint{N,CT},RT,Nv}
		return new{N,CP,RT,CT,Nv}(cpt)
	end
	function WaveletCollocationPoint{N,CP,RT,CT,Nv}(cpt::CP,children::SVector{N,SVector{2,WaveletCollocationPoint{N,CP,RT,CT,Nv}}},parents::SVector{N,WaveletCollocationPoint{N,CP,RT,CT,Nv}}) where {N,CT,CP<:AbstractCollocationPoint{N,CT},RT,Nv}
		return new{N,CP,RT,CT,Nv}(cpt,children,parents)
	end
end
end

@generated function invalid_collocation_point(::Type{HCP}) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}}
	icpt = HCP(invalid_collocation_point(CP))
	return quote
		$icpt
	end
end

@generated function root_point(::Type{HCP}) where {N,CP,RT,HCP<:AbstractHierarchicalCollocationPoint{N,CP,RT}}
	rp = root_point(CP)
	children = _unroll_(SVector{N,SVector{2,HCP}},SVector{2,HCP}(invalid_collocation_point(HCP),invalid_collocation_point(HCP)))
	parents = _unroll_(SVector{N,HCP},invalid_collocation_point(HCP))
	return quote
		$HCP($rp,$children,$parents)
	end
end

function HierarchicalCollocationPoint(::Type{Val{DIM}}, prnt::HCP, crd::CT, cp_interv::Tuple{CT,CT}, ptidx::Int) where {DIM,N,CT,CP<:AbstractCollocationPoint{N,CT},HCP<:AbstractHierarchicalCollocationPoint{N,CP}}
	return HCP(
				CP(Val{DIM}, cpt(prnt), crd, cp_interv, ptidx),
				_unroll_(SVector{N,SVector{2,HCP}},SVector{2,HCP}(invalid_collocation_point(HCP),invalid_collocation_point(HCP))),
				_unroll_(SVector{N,HCP},invalid_collocation_point(HCP),prnt,Val{DIM})
				)
end

include("./CollocationPoints/utils_HierarchicalCollocationPoints.jl")
include("./CollocationPoints/utils_WaveletCollocationPoints.jl")


