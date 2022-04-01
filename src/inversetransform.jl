using Distributions

function inversetransform(d::D, point_coord::CT) where {D<:Distribution{Univariate,Continuous}, CT <: Number}
	f = let d=d; x->quantile(d, x); end
	return inversetransform(f, point_coord)
end

function inverse_chebyshev_cdf(x)
	return -cos(π * x)
end

function inversetransform(invcdf::F, point_coord::CT) where {F<:Function, CT<:Number}
	sgn = sign(point_coord)
	absx = (abs(point_coord)+1)/2
	@assert .5 <= absx <= 1.
	return sgn * invcdf(absx)
end
