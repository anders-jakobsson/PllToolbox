struct PllStepResponse{T}
	time::Vector{T}
	response::Vector{T}
end

struct PllStepResponseWrapper{T} <: Tables.AbstractColumns
	sr::PllStepResponse{T}
end

Tables.istable(x::PllStepResponse) = true
Tables.columnaccess(::Type{PllStepResponse}) = true
Tables.columns(x::PllStepResponse) = PllStepResponseWrapper(x)
function Tables.columnnames(::PllStepResponseWrapper)
	names = uppercasefirst.(String.(fieldnames(PllStepResponse)))
	Symbol.(names)
end
function Tables.getcolumn(x::PllStepResponseWrapper, k::Int)
	sr = getfield(x, :sr)
	getfield(sr,k)
end
function Tables.getcolumn(x::PllStepResponseWrapper, name::Symbol)
	sr = getfield(x, :sr)
	getfield(sr,Symbol(lowercase(String(name))))
end
function Tables.schema(::PllStepResponseWrapper{T}) where T
	Tables.Schema(Tables.columnnames,(T,T))
end

@userplot StepPlot

stepplot(pll::PLL, tfinal=0, npoints=100; kw...) = stepplot(pllstep(pll,tfinal,npoints); kw...)
stepplot!(pll::PLL, tfinal=0, npoints=100; kw...) = stepplot!(pllstep(pll,tfinal,npoints); kw...)

@recipe function f(x::StepPlot)
	arg = x.args[1]
	t    = arg.time
	y    = arg.response
	tₘᵢₙ = minimum(t)
	tₘₐₓ = maximum(t)
	
	framestyle := :semi
	xlims      := (tₘᵢₙ,tₘₐₓ)
	seriestype := :line
	linewidth  --> 2
	legend     --> false
	grid       --> true
	
	@series begin
		xlabel  := "Time (s)"
		ylabel  := "Amplitude (radians or Hz)"
		t, y
	end
end




export pllstep

"""
	y,t = pllstep(pll:PLL, tfinal=0, npoints=100)

Return the step response of a PLL. Computes the step response at npoints between 0 and 
tfinal. If tfinal is non-positive, it will be calculated based on system dynamics.
"""
function pllstep(pll::PLL, tfinal=0, npoints=100)
	CL = pll["Closed-loop"]
	num = coeffs(CL.num)
	den = coeffs(CL.den)
	m = length(num)
	n = length(den)
	num /= den[end]
	den /= den[end]
	pnum = Polynomial(num)
	pden = Polynomial(den)
	if isnan(pnum(0)/pden(0))
		error("PLL is unstable, unable to compute step response.")
	end
	if tfinal<=0
		τ = minimum(abs.(real.(roots(pden))))
		tfinal = -log(1e-3)/τ
	end
	
	A = [zeros(n-2,1) I(n-2);-den[1:end-1]']
	B = [zeros(n-2,1);1]
	C = [num' zeros(1,n-m-1)]
	tspan = (0,tfinal)
	func(y,p,t) = A*y+B
	y0 = zeros(n-1,1)
	prob = ODEProblem(func,y0,tspan)
	sol = solve(prob,Rodas4P(), abstol=1e-12, reltol=1e-6)
	t = LinRange(0,tfinal,npoints)
	y = [(C*sol(tₖ))[1] for tₖ in t]
	return PllStepResponse(collect(t), y)
end