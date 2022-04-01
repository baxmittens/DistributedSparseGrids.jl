
function generate_unroll_code_dim(::Type{Val{N}}, ::Type{Val{DIM}}) where {N,DIM}
	insert_string(str::String,sub_str::String,i::Int,N::Int) = begin; str *= sub_str; i < N ? str*="," : nothing; return str; end
	insert_parentstr(str::String,i::Int,N::Int) = insert_string(str,"prnt[$i]",i,N)
	insert_elstr(str::String,i::Int,N::Int) = insert_string(str,"el",i,N)
	insert_next(str::String,i::Int,N::Int,DIM::Int) = i == DIM ? insert_elstr(str,i,N) : insert_parentstr(str,i,N)
	str = "("
	for i = 1:N
		str = insert_next(str,i,N,DIM)
	end
	str *= ")"
	return Meta.parse(str)
end

function generate_unroll_code_dim_const_el(::Type{Val{N}}, ::Type{Val{DIM}}) where {N,DIM}
	insert_string(str::String,sub_str::String,i::Int,N::Int) = begin; str *= sub_str; i < N ? str*="," : nothing; return str; end
	insert_parentstr(str::String,i::Int,N::Int) = insert_string(str,"const_el",i,N)
	insert_elstr(str::String,i::Int,N::Int) = insert_string(str,"el",i,N)
	insert_next(str::String,i::Int,N::Int,DIM::Int) = i == DIM ? insert_elstr(str,i,N) : insert_parentstr(str,i,N)
	str = "("
	for i = 1:N
		str = insert_next(str,i,N,DIM)
	end
	str *= ")"
	return Meta.parse(str)
end

function generate_unroll_code_const_el(::Type{Val{N}}) where {N}
	insert_string(str::String,sub_str::String,i::Int,N::Int) = begin; str *= sub_str; i < N ? str*="," : nothing; return str; end
	insert_parentstr(str::String,i::Int,N::Int) = insert_string(str,"const_el",i,N)
	insert_next(str::String,i::Int,N::Int) = insert_parentstr(str,i,N)
	str = "("
	for i = 1:N
		str = insert_next(str,i,N)
	end
	str *= ")"
	return Meta.parse(str)
end

function gen_unroll_dyn_code(::Type{Val{N}}) where {N}
	str = "if dim == 1
	return _unroll_(prnt,el,Val{1})"
	for i = 2:N
		str*="
elseif dim == $i
	return _unroll_(prnt,el,Val{$i})"
	end
	str *= "
else
	error()
end"
	return Meta.parse(str)
end

@generated function _unroll_(prnt::SVector{N,T}, el::T, ::Type{Val{DIM}}) where {N,T,DIM}
	code = generate_unroll_code_dim(Val{N}, Val{DIM})
	return quote
		SVector{$N,$T}($code)
	end
end

@generated function _unroll_dyn_dim_(prnt::SVector{N,T}, el::T, dim::Int) where {N,T}
	code = gen_unroll_dyn_code(Val{N})
	return quote
		$code
	end
end

@generated function _unroll_(::Type{SVector{N,T}}, const_el::T, el::T, ::Type{Val{DIM}}) where {N,T,DIM}
	code = generate_unroll_code_dim_const_el(Val{N}, Val{DIM})
	return quote
		SVector{$N,$T}($code)
	end
end
@generated function _unroll_(::Type{SVector{N,T}}, const_el::T) where {N,T}
	code = generate_unroll_code_const_el(Val{N})
	return quote
		SVector{$N,$T}($code)
	end
end