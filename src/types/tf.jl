using Polynomials

export Tf

"""
		Tf(num,den)
	
Create a TF transfer function object given numerator coefficients in num and denominator
coefficients in den. Both arguments must be iterable. 
"""
struct Tf
	num::Polynomial
	den::Polynomial

	function Tf(pnum::Polynomial,pden::Polynomial)
		pgcd = gcd(pnum,pden)
		return new(pnum÷pgcd,pden÷pgcd)
	end
end

function Tf(num,den)
	pnum = Polynomial(Iterators.reverse(num),:s)
	pden = Polynomial(Iterators.reverse(den),:s)
	return Tf(pnum,pden)
end


# Custom pretty-printing.
function Base.show(io::IO, x::Tf)
	printfun(io,x) = printpoly(io, x, descending_powers=true, mulsymbol="")
	numstr = sprint(printfun, x.num)
	denstr = sprint(printfun, x.den)
	nlen = length(numstr)
	dlen = length(denstr)
	npad = abs(nlen-dlen)÷2 + min(nlen,dlen)
	line = "-"^max(nlen,dlen)
	println(lpad(numstr,npad))
	println(line)
	println(lpad(denstr,npad))
	return nothing
end



function Base.:+(a::Tf, b::Tf)
	num = a.num*b.den + a.den*b.num
	den = a.den * b.den
	return Tf(num,den)
end


function Base.:-(a::Tf, b::Tf)
	num = a.num*b.den - a.den*b.num
	den = a.den * b.den
	return Tf(num,den)
end


function Base.:*(a::Tf, b::Tf)
	num = a.num * b.num
	den = a.den * b.den
	return Tf(num,den)
end


function Base.:/(a::Tf, b::Tf)
	num = a.num * b.den
	den = a.den * b.num
	return Tf(num,den)
end


function evalfr(tf::Tf, s::Complex)
	return tf.num(s) / tf.den(s)
end

function freqresp(tf::Tf, ω)
	return evalfr.(tf, 1im*ω)
end