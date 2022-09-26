
function num2si(x::Real)
	xs,xp,xk = num2si([x])
	return xs[1],xp,xk
end

function num2si(x::AbstractArray{<:Real})
	prf = ["y","z","a","f","p","n","μ","m","","k","M","G","T","P","E","Z","Y"]
	pow = [-24,-21,-18,-15,-12,-9,-6,-3,0,3,6,9,12,15,18,21,24]
	xs = copy(x)

	a = log10(minimum(abs.(x)))
	b = log10(maximum(abs.(x)))
	if any(isnan.(x)) || any(isinf.(x)) || (0<=a<3 && 0<=b<3)
		return xs,"",1.0
	end
	
	
	e = ((b-a)-3)/2
	o = -(a+e)
	imin = abs.(pow.+o) .== minimum(abs.(pow.+o))
	mi = findlast(imin)
	
	xk = 10.0^(-pow[mi])
	xp = prf[mi]
	xs = xs.*xk
	return xs,xp,xk
end

function num2si(::Type{String}, x::Integer)
	xs,xp, = num2si(x)
	space = isempty(xp) ? "" : " "
	return string(xs) * space * xp
end

function num2si(::Type{String}, x::Real, n::Integer=2)
	xs,xp, = num2si(x)
	space = isempty(xp) ? "" : " "
	return string(round(xs,digits=n)) * space * xp
end