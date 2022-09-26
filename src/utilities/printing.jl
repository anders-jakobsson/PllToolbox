function pllpolyrepr(poly::AbstractVector{T}) where {T<:Real}
	str = ""
	first = true
	len = length(poly)
	power = len-1
	for k in 1:len
		if poly[k]==0
			power -= 1
			continue
		end
		str *= first ? repr(poly[k]) : repr(abs(poly[k]))
		first = false
		if power>0
			str *= "s"
		end
		if power>1
			str *= "^$(power)"
		end
		if k<len && any(x->x!=0,poly[k+1:end])
			str *= sign(poly[k+1])>0 ? " + " : " - "
		end
		power -= 1
	end
	return str
end




# function plltfrepr(H::LTISystem)
# 	numstr = pllpolyrepr(num(H)[1])
# 	denstr = pllpolyrepr(den(H)[1])
# 	len = max(length(numstr), length(denstr))
# 	numstr = " "^max(0,round(Int,(len-length(numstr))/2)) * numstr
# 	denstr = " "^max(0,round(Int,(len-length(denstr))/2)) * denstr
# 	str = numstr * "\n" * "-"^len * "\n" * denstr
# 	return str
# end