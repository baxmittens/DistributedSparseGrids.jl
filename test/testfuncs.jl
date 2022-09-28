#From The Integration Test Problems Toolbox, Michael Baudin
function PRODONES(x::AV1,x0::AV2) where {CT1,CT2,AV1<:AbstractVector{CT1},AV2<:AbstractVector{CT2}}
	g = (x,x0)-> x<x0 ? -1.0 : 1.0
	return foldl(*,map(g,x,x0))
end

function NCUBE(x::AV1) where {CT1,AV1<:AbstractVector{CT1}}
	g = x-> x < -0.5+pi/10 ? 0.0 : (x <= 0.5+pi/10 ? 1.0 : 0.0)
	return foldl(*,map(g,x))
end

function NBALL(x::AV1) where {CT1,AV1<:AbstractVector{CT1}}
	return sqrt(sum(x.^2))<1.0
end

function PRODEXP(x::AV1,x0::AV2)  where {CT1,CT2,AV1<:AbstractVector{CT1},AV2<:AbstractVector{CT2}}
	xx = (x.+1.0)./2.0
	N = length(x)
	w = ((15.0*exp(15.0)+15)/(13.0*exp(15.0)+17))^0.5
	g = (_x,_x0) -> (exp(30*(_x-_x0)-15)-1)/(exp(30*(_x-_x0)-15)+1)
	return w^N*foldl(*,map(g,xx,x0))
end

function HELLEKALEK(x::AV1,x0::AV2,α::Int)  where {CT1,CT2,AV1<:AbstractVector{CT1},AV2<:AbstractVector{CT2}}
	xx = (x.+1.0)./2.0
	γ = sqrt(α^2/((2*α+1))*(α+1)^2)
	h = x -> x^α - 1.0/(α+1)
	g = x -> h(x)/γ
	return foldl(*,map(g,xx-x0))
end

function ROOSARNOLD1(x::AV1,x0::AV2)  where {CT1,CT2,AV1<:AbstractVector{CT1},AV2<:AbstractVector{CT2}}
	xx = (x.+1.0)./2.0
	N = length(x)
	v = 1/(3*N)
	g = x -> abs(4*x-2)/N
	return 1/sqrt(v)*sum(map(g,xx-x0)-1.0)
end

function ROOSARNOLD2(x::AV1,x0::AV2)  where {CT1,CT2,AV1<:AbstractVector{CT1},AV2<:AbstractVector{CT2}}
	xx = (x.+1.0)./2.0
	N = length(x)
	v = (4.0/3.0)^N-1.0
	g = _x -> abs(4*_x-2)
	return 1/sqrt(v)*foldl(*,map(g,xx-x0)-1.0)
end

function RST(x::AV1,x0::AV2,αs::AV3)  where {CT1,CT2,AV1<:AbstractVector{CT1},AV2<:AbstractVector{CT2},AV3<:AbstractVector{Int}}
	xx = (x.+1.0)./2.0
	g = (_x,_a) -> (abs(4*_x-2)+_a)/(1.0+_a)
	γ = _a -> sqrt(1.0/(3.0*(1+_a)^2))
	v = foldl(*,1.0.+map(γ,αs).^2)-1.0
	return 1/sqrt(v)*foldl(*,map(g,xx.-x0,αs).-1.0)
end

function RST1(x::AV1,x0::AV2) where {CT1,CT2,AV1<:AbstractVector{CT1},AV2<:AbstractVector{CT2}}
	a = [1 for i = 1:length(x)]
	return RST(x,x0,a)
end

function RST2(x::AV1,x0::AV2) where {CT1,CT2,AV1<:AbstractVector{CT1},AV2<:AbstractVector{CT2}}
	a = [i for i = 1:length(x)]
	return RST(x,x0,a)
end

function RST3(x::AV1,x0::AV2) where {CT1,CT2,AV1<:AbstractVector{CT1},AV2<:AbstractVector{CT2}}
	a = [i^2 for i = 1:length(x)]
	return RST(x,x0,a)
end

function MC_integral(::Type{Val{N}},::Type{Val{S}},f::F) where {N,S,F<:Function}
	res = BigFloat(0.0)
	for i = 1:S
		x = rand(BigFloat, N).*2.0.-1.0
		res += f(x)
	end
	return 4.0*res/S
end

function testfunc(x::AV1) where {CT1,AV1<:AbstractVector{CT1}}
	return 1.0/(abs(foldl(-, (x.-1.0).^2, init=sqrt(2.0)))+0.5)
	#return 1.0/(abs(sqrt(2.0)-x[1]^2-x[2]^2)+0.5)
end


