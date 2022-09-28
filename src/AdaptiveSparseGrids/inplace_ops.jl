import LinearAlgebra

add!(a::Matrix{Float64}, b::Matrix{Float64}) = a .+= b
add!(a::Matrix{Float64}, b::Float64) = a .+= b
LinearAlgebra.mul!(a::Matrix{Float64}, b::Float64) = a .*= b
LinearAlgebra.mul!(a::Matrix{Float64}, b::Matrix{Float64}, c::Float64) = begin; for i = 1:length(a); a[i] = b[i]*c; end; return nothing; end
minus!(a::Matrix{Float64}, b::Matrix{Float64}) = a .-= b
minus!(a::Matrix{Float64}, b::Float64) = a .-= b
pow!(a::Matrix{Float64}, b::Int64) = a .^= b
max!(a::Matrix{Float64}, b::Matrix{Float64}) = a .= max.(a,b)
min!(a::Matrix{Float64}, b::Matrix{Float64}) = a .= min.(a,b)