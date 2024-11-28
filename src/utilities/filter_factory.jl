function pllfilter(name::String, k::Number, z::Vector{<:Number}, p::Vector{<:Number})
	Block(name, zpk(z,p,k), AbstractNoise[])
end








""""
    pllfilter(name::String, Aloop::Real, fc::Real, pm::Real, R₂₁::Real=0, R₃₂::Real=0)

Construct a passive PLL filter block from the given corner frequency and pole ratios etc.


"""
function pllfilter(name::String, Aloop::Real, fc::Real, pm::Real, temp::Real, r₂₁::Real=0, r₃₂::Real=0)
	ωc = 2π*fc
	ϕm = π*pm/180
	if iszero(r₂₁)
		R,C = _pllfilter2(Aloop, ωc, ϕm)
	elseif iszero(r₃₂)
		~,R,C = _pllfilter3(Aloop, ωc, ϕm, r₂₁, r₃₂)
	else
		R,C = _pllfilter4(Aloop, ωc, ϕm, r₂₁, r₃₂)
	end
	pllfilter(name, R, C, temp)
end


function _pllfilter2(Aloop, ωc, ϕm)
	p₁ = (sqrt(tan(ϕm)^2+1)-tan(ϕm)) / ωc
	z  = 1/(ωc^2*p₁)
	A₀ = (Aloop/ωc^2) * sqrt((1+(ωc*z)^2) / (1+(ωc*p₁)^2))
	C₁ = A₀*p₁/z
	C₂ = A₀-C₁
	R₂ = z/C₂
	return [R₂,],[C₁,C₂]
end


function _pllfilter3(Aloop, ωc, ϕm, r₂₁, r₃₂)
	# The transfer function phase is given by:
	#    arg(H(jω)) = ϕ(ω) = atan(ω*z)-atan(ω*p₁)-atan(ω*p₁*r₂₁)-atan(ω*p₁*r₂₁*r₃₂)
	#
	# Finding the pole(s) and the zero involves maximizing the phase margin at ωc, 
	# that is:
	#    ϕ(ωc)-ϕm = 0
	#    dϕ(ωc)/dω = 0
	#
	# These two equations can be used to numerically solve for z and p₁.
	function f!(F,x)
		zf,pf = x
		F[1] = atan(ωc*zf)-atan(ωc*pf)-atan(ωc*pf*r₂₁)-atan(ωc*pf*r₂₁*r₃₂)-ϕm
		F[2] = zf/(1+(ωc*zf)^2)-pf/(1+(ωc*pf)^2)-pf*r₂₁/(1+(ωc*pf*r₂₁)^2)-pf*r₂₁*r₃₂/(1+(ωc*pf*r₂₁*r₃₂)^2)
	end

	pinit = (sec(ϕm)-tan(ϕm)) / (ωc*(1+r₂₁))
	zinit = 1/(ωc^2*pinit*(1+r₂₁))
	solution = nlsolve(f!, [zinit,pinit], ftol=1e-9, autodiff=:forward)
	
	z,p₁ = solution.zero
	p₂ = p₁*r₂₁

	A₀ = (Aloop/ωc^2) * sqrt((1+(ωc*z)^2) / ((1+(ωc*p₁)^2)*(1+(ωc*p₂)^2)))
	A₁ = A₀*(p₁+p₂)
	A₂ = A₀*p₁*p₂
	C₁ = (A₂/z^2)*(1+sqrt(1+(z/A₂)*(z*A₀-A₁)))
	C₃ = (z*A₁*C₁-(z*C₁)^2-A₀*A₂) / (C₁*z^2-A₂)
	C₂ = A₀-C₁-C₃
	R₂ = z/C₂
	R₃ = A₂/(C₁*C₃*z)

	return (z,p₁,p₂,A₀),[R₂,R₃],[C₁,C₂,C₃]
end


function _pllfilter4(Aloop, ωc, ϕm, r₂₁, r₃₂)
	# Solution for a 3rd-order filter is used as a starting point
	sol,RS,CS = _pllfilter3(Aloop,ωc,ϕm,r₂₁,r₃₂)
	(z,p₁,p₂,Aₛ) = sol
	p₃ = p₁*r₂₁*r₃₂
	A̅₁ = Aₛ*(p₁+p₃)
	A̅₂ = Aₛ*p₁*p₃
	C̅₁ = (A̅₂/z^2)*(1+sqrt(1+(z/A̅₂)*(z*Aₛ-A̅₁)))
	C̅₃ = (z*A̅₁*C̅₁-(z*C̅₁)^2-Aₛ*A̅₂) / (C̅₁*z^2-A̅₂)
	R̅₃ = A̅₂/(C̅₁*C̅₃*z)
	C₁ = (C₁+C̅₁)/2
	R₃ = (R₃+R̅₃)/2
	A₀ = (Kl/ωc^2) * sqrt((1+(ωc*z)^2) / ((1+(ωc*p₁)^2)*(1+(ωc*p₂)^2)*(1+(ωc*p₃)^2)))
	A₁ = A₀*(p₁+p₂+p₃)
	A₂ = A₀*(p₁*p₂ + p₁*p₃ + p₂*p₃)
	A₃ = A₀*p₁*p₂*p₃
	k₀ = A₂/A₃ - 1/z - 1/(C₁*R₃) - (A₀-C₁)*z*R₃*C₁/A₃
	k₁ = A₁ - z*A₀ - A₃/(z*R₃*C₁) - (A₀-C₁)*R₃*C₁
	a = A₃/(z*C₁)^2
	b = z + R₃*(C₁-A₀) + A₃*(1/z-k₀)/(z*C₁)
	c = k₁ - k₀*A₃/z
	C₂ = (sqrt(b^2-4*a*c)-b)/(2*a)
	C₃ = z*A₃*C₁/(R₃*(k₀*z*A₃*C₁ - C₂*(A₃-R₃*(z*C₁)^2)))
	C₄ = A₀ - C₁ - C₂ - C₃
	R₂ = z/C₂
	R₄ = A₃/(z*R₃*C₁*C₃*C₄)

	return [R₂,R₃,R₄],[C₁,C₂,C₃,C₄]
end







function pllfilter(name::String, res::AbstractVector{<:Real}, cap::AbstractVector{<:Real}, temp::Real=300)
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

	Block(name, H, noise),R,C
end
