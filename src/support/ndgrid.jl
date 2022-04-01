function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end

function ndgrid(vs::AbstractVector{T}...) where {T}
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array{T}(undef, sz), n)
    s = 1
    for i=1:n
        a = out[i]::Array{T}
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end
