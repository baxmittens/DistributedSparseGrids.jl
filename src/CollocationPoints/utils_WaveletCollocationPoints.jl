function set_wavelet_coeffs!(wcpt::WaveletCollocationPoint{N,CP,RT,CT,Nv},wcoeffs::AV1) where {N,CT,CP,RT,Nv,AV1<:AbstractVector{CT}}
	@assert Nv >= length(wcoeffs) "Number of wavelet coefficients is $(length(wcoeffs)) but should be at most $(Nv)."
	tr_zeros = zeros(CT, Nv-length(wcoeffs))
	wcpt.wavelet_coeffs = SVector{Nv,CT}(wcoeffs...,tr_zeros...)
	return nothing
end

wavelet_coeffs(wcpt::WaveletCollocationPoint) = wcpt.wavelet_coeffs

function Base.isequal(x::WaveletCollocationPoint,y::WaveletCollocationPoint)
	return all(i_multi(x).==i_multi(y)) && all(pt_idx(x).==pt_idx(y))
end

