using Polynomials


"""
	TF
	
Rational linear, time invariant (LTI) transfer function data type. 
"""
struct TF
	num::Polynomial
	den::Polynomial
	name::String
	type::Symbol
end


"""
	TF(pnum::Polynomial{T,:s}, pden::Polynomial{T,:s}, name::String="")

Create a transfer function with the given name and numerator/denominator polynomial.
"""
function TF(pnum::Polynomial, pden::Polynomial, name::String="")
	num = coeffs(pnum)
	den = coeffs(pden)
	npr = min(findfirst(!iszero, num), findfirst(!iszero, den))
	pnm = Polynomial(num[npr:end],:s)
	pdn = Polynomial(den[npr:end],:s)
	TF(pnm, pdn, name, classify(pnm,pdn))
end


"""
	TF(num=1, den=1, name::String="")

Create a transfer function with the given name and numerator/denominator coefficients.
Arguments `num` and `den` are the coefficients in descending orders of `s`, and must be 
a numeric array, tuple or scalar.
"""
function TF(num=1, den=1, name::String="")
	pnum = Polynomial(Iterators.reverse(_coeff_convert(num)),:s)
	pden = Polynomial(Iterators.reverse(_coeff_convert(den)),:s)
	return TF(pnum, pden, name) 
end


"""
    TF(orig::TF, num=x.num, den=orig.den, name=orig.name)

Create a copy of transfer function `orig`, optionally overiding the fields with new values.
"""
TF(orig::TF; num=reverse(coeffs(orig.num)), den=reverse(coeffs(orig.den)), name=orig.name) = TF(num,den,name)


# Coefficient connversion:
_coeff_convert(x::Number) = [x]
_coeff_convert(x::Tuple) = _coeff_convert(collect(x))
_coeff_convert(x::AbstractArray{T,N}) where {T<:Number,N} = x




# Custom pretty-printing:
Base.show(io::IO, x::TF) = print(io, "$(x.name) ($(x.type))")
Base.show(io::IO, ::MIME"text/plain",  x::TF) = _print_tf(io, x, MIME("text/plain"),  "-")
Base.show(io::IO, ::MIME"text/html",  x::TF)  = _print_tf(io, x, MIME("text/html"),   "─")
function _print_tf(io::IO, x::TF, mimetype, linechar="-")
	printfun(io,x) = printpoly(io, x, mimetype, descending_powers=true, compact=true, mulsymbol="")
	numstr = sprint(printfun, x.num)
	denstr = sprint(printfun, x.den)
	nlen = length(numstr)
	dlen = length(denstr)
	line = linechar^max(nlen,dlen)
	name = isempty(x.name) ? "" : x.name*": "
	npad = length(name) + max(dlen-nlen,0)÷2 + nlen
	dpad = length(name) + max(nlen-dlen,0)÷2 + dlen
	println(io, lpad(numstr,npad))
	println(io, name*line)
	print(io, lpad(denstr,dpad))
	return nothing
end
Base.show(io::IO, ::MIME"text/latex",  x::TF)  = println(io, raw"\frac{"*printpoly(io, x.num, MIME("text/latex"))*"}{"*printpoly(io, x.den, MIME("text/latex"))*"}")


# Conversion and promotion:
Base.convert(::Type{TF}, x::Number) = TF(x)
Base.promote_rule(::Type{TF}, ::Type{T}) where {T<:Number} = TF
Base.:+(a::Union{T1,T2},b::Union{T1,T2}) where {T1<:Number,T2<:TF} = +(promote(a,b)...)
Base.:-(a::Union{T1,T2},b::Union{T1,T2}) where {T1<:Number,T2<:TF} = -(promote(a,b)...)
Base.:*(a::Union{T1,T2},b::Union{T1,T2}) where {T1<:Number,T2<:TF} = *(promote(a,b)...)
Base.:/(a::Union{T1,T2},b::Union{T1,T2}) where {T1<:Number,T2<:TF} = /(promote(a,b)...)




# Arithmetic:
function Base.:-(a::TF)
	return TF(-a.num,a.den,a.name)
end


function Base.:+(a::TF, b::TF)
	num = a.num*b.den + a.den*b.num
	den = a.den * b.den
	return TF(num,den,"")
end


function Base.:-(a::TF, b::TF)
	num = a.num*b.den - a.den*b.num
	den = a.den * b.den
	return TF(num,den,"")
end


function Base.:*(a::TF, b::TF)
	num = a.num * b.num
	den = a.den * b.den
	return TF(num,den,"")
end


function Base.:/(a::TF, b::TF)
	num = a.num * b.den
	den = a.den * b.num
	return TF(num,den,"")
end





"""
    (tf::TF)(s::Complex)

Evaluate the transfer function tf at s.
"""
function (tf::TF)(s::Complex)
	return tf.num(s) / tf.den(s)
end




"""
	m,phi = bode(tf::TF, f::AbstractArray{<:Real,N}; magunit=:dB, phaseunit=:rad)

Return the magnitude response `m` and phase response `phi` of transfer function `tf` at the
frequencies in `f` (in Hertz). The units are set by `magunit` (:lin/:dB) and `phaseunit` 
(:deg/:rad). 
"""
function bode(tf::TF, f::AbstractArray{T,N}; magunit::Symbol=:dB, phaseunit::Symbol=:deg) where {T<:Real, N}
	magunit ∈ (:lin,:dB) || error("Unknown magnitude unit `$(magunit)`")
	phaseunit ∈ (:rad,:deg) || error("Unknown phase unit `$(phaseunit)`")
	s = 2im*π*f
	r = tf.num.(s) ./ tf.den.(s)
	m = magunit===:lin ? abs.(r) : dB20.(abs.(r))
	ϕ = phaseunit===:rad ? unwrap(angle.(r)) : rad2deg.(unwrap(angle.(r)))
	(m,ϕ)
end




"""
	m = magresp(tf::TF, f::AbstractArray{<:Real,N})

Return the magnitude response of transfer function `tf` at frequencies in `f` (in Hertz).
"""
function magresp(tf::TF, f::AbstractArray{T,N}) where {T<:Real, N}
	s = 2im*π*f
	r = tf.num.(s) ./ tf.den.(s)
	abs.(r)
end




"""
	phi = phaseresp(tf::TF, f::AbstractArray{<:Real,N})

Return the phase response in radians of transfer function `tf` at frequencies `f` (in Hertz).
"""
function phaseresp(tf::TF, f::AbstractArray{T,N}) where {T<:Real, N}
	s = 2im*π*f
	r = tf.num.(s) ./ tf.den.(s)
	unwrap(angle.(r))
end




"""
    unwrap(ϕ, period=2π)

Unwrap the phase `ϕ∈[0,2π]`, removing instantaneous phase jump.
For the in-place version of this function, see [`unwrap!`](@ref).
"""
function unwrap(ϕ, period=2π)
	p = convert(eltype(ϕ), period)
	v = first(ϕ)
	Φ = similar(ϕ)
	@inbounds for k = eachindex(ϕ)
		Φ[k] = v = v + rem(ϕ[k]-v, p, RoundNearest)
	end
	return Φ
end


"""
    unwrap!(ϕ, period=2π)

Unwrap the phase `ϕ∈[0,2π]`, removing instantaneous phase jump.
"""
function unwrap!(ϕ, period=2π)
	p = convert(eltype(ϕ), period)
	v = first(ϕ)
	@inbounds for k = eachindex(ϕ)
		ϕ[k] = v = v + rem(ϕ[k]-v, p, RoundNearest)
	end
end



export classify
classify(tf::TF) = classify(tf.num,tf.den)
function classify(pnum,pden)
	nₘₐₓ = length(pnum)-1
	nₘᵢₙ = findfirst(!iszero, pnum)
	dₘₐₓ = length(pden)-1
	dₘᵢₙ = findfirst(!iszero, pden)
	
	nₘₐₓ==dₘₐₓ==0 && return :constant
	nₘₐₓ==dₘₐₓ && return :highpass
	nₘᵢₙ>0 && return :bandpass
	return :lowpass
end