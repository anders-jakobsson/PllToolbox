export pllfilter



function _rms(x::Vector{Float64})
	sqrt(sum(x.^2)/length(x))
end




function pllfilter(name::String, k::Number, z::Vector{<:Number}, p::Vector{<:Number})
	Block(name, zpk(z,p,k), AbstractNoise[])
end








""""
    pllfilter(name::String, Aloop::Real, fc::Real, pm::Real, R₂₁::Real=0, R₃₂::Real=0)

Construct a passive PLL filter block from the given corner frequency and pole ratios etc.


"""
function pllfilter(name::String, Aloop::Real, fc::Real, pm::Real, temp::Real, R₂₁::Real=0, R₃₂::Real=0)
	Kl  = Float64(Aloop)
	ωc  = 2*pi*Float64(fc)
	ϕm  = pi*Float64(pm)/180
	r21 = Float64(R₂₁)
	r32 = Float64(R₃₂)

	# Second order filter has a closed form. Higher order filters are found numerically.
	if r21==0
		p1 = (sqrt(tan(ϕm)^2+1)-tan(ϕm)) / ωc
		z  = 1/(ωc^2*p1)
		A0 = (Kl/ωc^2) * sqrt((1+(ωc*z)^2) / (1+(ωc*p1)^2))
		C₁ = A0*p1/z
		C₂ = A0-C₁
		C₃ = C₄ = 0
		R₂ = z/C₂
		R₃ = R₄ = 0

	else
		# The transfer function phase is given by:
		#    arg(H(jω)) = ϕ(ω) = atan(ω*z)-atan(ω*p1)-atan(ω*p1*r21)-atan(ω*p1*r21*r32)
		#
		# Finding the pole(s) and the zero involves maximizing the phase margin at ωc, 
		# that is:
		#    ϕ(ωc)-ϕm = 0
		#    dϕ(ωc)/dω = 0
		#
		# These two equations can be used to numerically solve for z and p1. This is done 
		# by forming a vector valued function and its Jacobian matrix.

		ϕc = (z,p1) -> atan(ωc*z)-atan(ωc*p1)-atan(ωc*p1*r21)-atan(ωc*p1*r21*r32)-ϕm
		Dωϕc = (z,p1) -> z/(1+(ωc*z)^2)-p1/(1+(ωc*p1)^2)-p1*r21/(1+(ωc*p1*r21)^2)-p1*r21*r32/(1+(ωc*p1*r21*r32)^2)
		f = x::Vector{Float64} -> [
				ϕc(x[1],x[2])
				Dωϕc(x[1],x[2])
			]

		J11 = x::Vector{Float64} -> ωc/(1+(ωc*x[1])^2)
		J12 = x::Vector{Float64} -> -ωc/(1+(ωc*x[2])^2) - ωc*r21/(1+(ωc*r21*x[2])^2) - ωc*r21*r32/(1+(ωc*r21*r32*x[2])^2)
		J21 = x::Vector{Float64} -> (1-(ωc*x[1])^2) / (1+(ωc*x[1])^2)^2
		J22 = x::Vector{Float64} -> -(1-(ωc*x[2])^2)/(1+(ωc*x[2])^2)^2 - r21*(1-(ωc*r21*x[2])^2)/(1+(ωc*r21*x[2])^2)^2 - r21*r32*(1-(ωc*r21*r32*x[2])^2)/(1+(ωc*r21*r32*x[2])^2)^2
		J = x::Vector{Float64} -> [J11(x) J12(x);J21(x) J22(x)]

		p1init = (sec(ϕm)-tan(ϕm)) / (ωc*(1+r21))
		x = [1/(ωc^2*p1init*(1+r21)); p1init]
		for k=1:100
			x = x - J(x)\f(x)
			if _rms(f(x))<10*eps(Float64)
				break
			end
		end
		
		z = x[1]
		p1 = x[2]
		p2 = p1*r21
		p3 = p1*r21*r32

		A0 = (Kl/ωc^2) * sqrt((1+(ωc*z)^2) / ((1+(ωc*p1)^2)*(1+(ωc*p2)^2)))
		A1 = A0*(p1+p2)
		A2 = A0*p1*p2
		C₁ = (A2/z^2)*(1+sqrt(1+(z/A2)*(z*A0-A1)))
		C₃ = (z*A1*C₁-(z*C₁)^2-A0*A2) / (C₁*z^2-A2)
		C₂ = A0-C₁-C₃
		C₄ = 0
		R₂ = z/C₂
		R₃ = A2/(C₁*C₃*z)
		R₄ = 0
		if p3>0
			A̅1 = A0*(p1+p3)
			A̅2 = A0*p1*p3
			C̅₁ = (A̅2/z^2)*(1+sqrt(1+(z/A̅2)*(z*A0-A̅1)))
			C̅₃ = (z*A̅1*C̅₁-(z*C̅₁)^2-A0*A̅2) / (C̅₁*z^2-A̅2)
			R̅₃ = A̅2/(C̅₁*C̅₃*z)
			C₁ = (C₁+C̅₁)/2
			R₃ = (R₃+R̅₃)/2
			A0 = (Kl/ωc^2) * sqrt((1+(ωc*z)^2) / ((1+(ωc*p1)^2)*(1+(ωc*p2)^2)*(1+(ωc*p3)^2)))
			A1 = A0*(p1+p2+p3)
			A2 = A0*(p1*p2 + p1*p3 + p2*p3)
			A3 = A0*p1*p2*p3
			k0 = A2/A3 - 1/z - 1/(C₁*R₃) - (A0-C₁)*z*R₃*C₁/A3
			k1 = A1 - z*A0 - A3/(z*R₃*C₁) - (A0-C₁)*R₃*C₁
			a = A3/(z*C₁)^2
			b = z + R₃*(C₁-A0) + A3*(1/z-k0)/(z*C₁)
			c = k1 - k0*A3/z
			C₂ = (sqrt(b^2-4*a*c)-b)/(2*a)
			C₃ = z*A3*C₁/(R₃*(k0*z*A3*C₁ - C₂*(A3-R₃*(z*C₁)^2)))
			C₄ = A0 - C₁ - C₂ - C₃
			R₂ = z/C₂
			R₄ = A3/(z*R₃*C₁*C₃*C₄)
		end
	end
	R = [R₂,R₃,R₄]
	C = [C₁,C₂,C₃,C₄]
	pllfilter(name, R, C, temp)
end







function pllfilter(name::String, res::Vector{<:Real}, cap::Vector{<:Real}, temp::Real=300, topo::Int=1)
	R = [convert(Vector{Float64},res); zeros(Float64, 3-length(res))]
	C = [convert(Vector{Float64},cap); zeros(Float64, 4-length(cap))]
	kT4 = 1.380649e-23*Float64(temp)*4
	ir0 = findlast(R.!=0)
	ic0 = findlast(C.!=0)
	if isnothing(ir0)
		error("All resistor values are zero")
	end
	if isnothing(ic0)
		error("All capacitor values are zero")
	end
	order = min(ir0::Int+1, ic0::Int)
	
	if topo==1
		R₂,R₃,R₄ = R[1],R[2],R[3]
		C₁,C₂,C₃,C₄ = C[1],C[2],C[3],C[4]
		a = C₁*C₂*C₃*C₄*R₂*R₃*R₄
		b = C₁*C₂*R₂*R₃*(C₃+C₄) + C₄*R₄*(C₂*R₂*(C₁+C₃)+C₃*R₃*(C₁+C₂))
		c = C₂*R₂*(C₁+C₃+C₄) + R₃*((C₁+C₂)*(C₃+C₄)) + C₄*R₄*(C₁+C₂+C₃)
		d = C₁+C₂+C₃+C₄
		t = C₂*R₂
		u = C₁*C₂*R₂
		v = C₁+C₂
		w = C₁*C₂*C₃*R₂*R₃
		x = C₂*R₂*(C₁+C₃) + C₃*R₃*(C₁+C₂)
		y = C₁+C₂+C₃
		fden = [a,b,c,d]
		H = TF([t,1]/d, [fden;0]/d, name)
		Hn2 = TF(C₂, fden)
		Hn3 = TF([u,v], fden)
		Hn4 = TF([w,x,y], fden)
		noise = WhiteNoise[]
		if order>1
			noise = [WhiteNoise("R2 noise", kT4*R₂, Hn2)]
		end
		if order>2
			noise = [noise;WhiteNoise("R3 noise", kT4*R₃, Hn3)]
		end
		if order>3
			noise = [noise;WhiteNoise("R4 noise", kT4*R₄, Hn4)]
		end

	elseif topo==2
		R₁,R₂ = R[1],R[2]
		C₁,C₂ = C[1],C[2]
		a = C₁*C₂*R₁*R₂
		b = R₁*C₁ + R₂*C₂ + R₁*C₂
		t = C₂*R₂
		u = C₁*C₂*R₂
		v = C₁+C₂
		fden = [a,b,1]
		H = TF([R₁*t,R₁], fden, name)
		Hn1 = TF(1, fden)
		Hn2 = TF([u,v], fden)
		noise = AbstractNoise[]

	elseif topo==3
		R₁,R₂ = R[1],R[2]
		C₁,C₂ = C[1],C[2]
		a = C₁*C₂*R₁*R₂
		b = R₁*C₁ + R₂*C₂ + R₁*C₂
		t = C₂*R₂
		u = C₁*C₂*R₂
		v = C₁+C₂
		fden = [a,b,1]
		H = TF([R₁*t,R₁], fden), name
		Hn1 = TF(1, fden)
		Hn2 = TF([u,v], fden)
		noise = AbstractNoise[]


	else

	end

	Block(name, H, noise),R,C
end




