cpt(cpt::CollocationPoint) = cpt
i_multi(cpt::CollocationPoint) = cpt.i_multi
i_multi(cpt::CollocationPoint,dim::Int) = cpt.i_multi[dim]
pt_idx(cpt::CollocationPoint) = cpt.pt_idx
pt_idx(cpt::CollocationPoint,dim::Int) = cpt.pt_idx[dim]
coords(cpt::CollocationPoint) = cpt.coords
coord(cpt::CollocationPoint,dim::Int) = cpt.coords[dim]
intervals(cpt::CollocationPoint) = cpt.interv
interval(cpt::CollocationPoint,dim::Int) = cpt.interv[dim]

getkey(cpt::CP) where {CP<:AbstractCollocationPoint} = (pt_idx(cpt),i_multi(cpt))
level(cpt::CP) where {N,CT<:Real,CP<:AbstractCollocationPoint{N,CT}} = sum(i_multi(cpt))-N+1
level(cpt::CP, dim::Int) where {CP<:AbstractCollocationPoint} = i_multi(cpt,dim)
isroot(cpt::CP,dim::Int) where {CP<:AbstractCollocationPoint} = pt_idx(cpt,dim) == i_multi(cpt,dim) == 1

Base.length(cpt::CP) where {CP<:AbstractCollocationPoint} = 1
Base.isvalid(cpt::CP) where {CP<:AbstractCollocationPoint} = level(cpt) > 0
Base.hash(cpt::CP) where {CP<:AbstractCollocationPoint} = hash(getkey(cpt))
#Base.isvalid(::Nothing) = false
hasvalid(cpts::SVector{2,CP}) where {CP<:AbstractCollocationPoint} = level(cpts[1]) > 0 || level(cpts[2]) > 0

function Base.string(cpt::CollocationPoint{N,CT}) where {N,CT<:Real}
	str = "CollocationPoint{$N,$CT}(level=$(level(cpt)), id=$(getkey(cpt)), x=$(coords(cpt))"
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

const Cpt = CollocationPoint
const CptID{N} = Tuple{SVector{N,Int},SVector{N,Int}}
level(id::CptID{N}) where {N} = sum(id[2])-N+1