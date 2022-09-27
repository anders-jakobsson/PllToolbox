struct PllBodeInfo
	fϕ::Float64
	ϕₘ::Float64
	fg::Float64
	gₘ::Float64
end

function Base.show(io::IO, bi::PllBodeInfo)
	pmfstr = num2si(String,bi.fϕ,1)
	pmystr = num2si(String,bi.ϕₘ,1)
	gmfstr = num2si(String,bi.fg,1)
	gmystr = num2si(String,bi.gₘ,1)
	print(io, "Phase/gain margin: $(pmystr)° at $(pmfstr)Hz / $(gmystr)dB at $(gmfstr)Hz")
end

function Base.show(io::IO, ::MIME"text/plain", bi::PllBodeInfo)
	pmfstr = num2si(String,bi.fϕ,1)
	pmystr = num2si(String,bi.ϕₘ,1)
	gmfstr = num2si(String,bi.fg,1)
	gmystr = num2si(String,bi.gₘ,1)
	println(io, "         Phase margin: $(pmystr)°")
	println(io, " Unity-gain frequency: $(pmfstr)Hz")
	println(io, "          Gain margin: $(gmystr)dB")
	  print(io, "Unity-phase frequency: $(gmfstr)Hz")
end


struct PllBodeResponse{T}
	frequency::Vector{T}
	magnitude::Vector{T}
	phase::Vector{T}
	info::PllBodeInfo
end

function Base.show(io::IO, x::PllBodeResponse)
	len = length(x.frequency)
	pm  = num2si(String,x.info.ϕₘ,1)
	if len==0
		print(io, "Empty Bode response with a $(pm)° phase margin")
	elseif len==1
		print(io, "Scalar Bode response with a $(pm)° phase margin")
	else
		print(io, "$len-point Bode response with a $(pm)° phase margin")
	end
end

function Base.show(io::IO, ::MIME"text/plain", x::PllBodeResponse)
	len  = length(x.frequency)
	if len==0
		println(io, "Empty Bode response:")
	elseif len==1
		println(io, "Scalar Bode response:")
	else
		println(io, "$len-point Bode response:")
	end
	
	fmin = num2si(String,x.frequency[1],1)
	fmax = num2si(String,x.frequency[end],1)
	println(io, "      Frequency range: $(fmin)Hz to $(fmax)Hz")

	mmin = num2si(String,minimum(x.magnitude),1)
	mmax = num2si(String,maximum(x.magnitude),1)
	println(io, "      Magnitude range: $(mmin)dB to $(mmax)dB")

	pmin = num2si(String,minimum(x.phase),1)
	pmax = num2si(String,maximum(x.phase),1)
	println(io, "          Phase range: $(pmin)° to $(pmax)°")

	show(io, MIME("text/plain"), x.info)
end


struct PllBodeResponseWrapper{T} <: Tables.AbstractColumns
	br::PllBodeResponse{T}
end

Tables.istable(x::PllBodeResponse) = true
Tables.columnaccess(::PllBodeResponse) = true
Tables.columns(x::PllBodeResponse) = PllBodeResponseWrapper(x)
function Tables.columnnames(::PllBodeResponseWrapper)
	Symbol.(uppercasefirst.(String.(fieldnames(PllBodeResponse)[1:end-1])))
end
function Tables.getcolumn(x::PllBodeResponseWrapper, k::Int)
	br = getfield(x, :br)
	getfield(br, k)
end
function Tables.getcolumn(x::PllBodeResponseWrapper, name::Symbol)
	br = getfield(x, :br)
	getfield(br, Symbol(lowercase(String(name))))
end
function Tables.schema(x::PllBodeResponseWrapper{T}) where T
	Tables.Schema(Tables.columnnames(x),(T,T,T))
end




@userplot BodePlot

bodeplot(pll::PLL; kw...) = bodeplot(pll,_bode_fvec(pll); kw...)
bodeplot(pll::PLL, f::AbstractVector{T}; kw...) where T<:Number = bodeplot(pllbode(pll,f); kw...)
bodeplot!(pll::PLL; kw...) = bodeplot!(pll,_bode_fvec(pll); kw...)
bodeplot!(pll::PLL, f::AbstractVector{T}; kw...) where T<:Number = bodeplot!(pllbode(pll,f); kw...)

@recipe function f(x::BodePlot)
	arg = x.args[1]
	f    = arg.frequency
	m    = arg.magnitude
	ϕ    = arg.phase
	fₘᵢₙ = minimum(f)
	fₘₐₓ = maximum(f)
	fₘ   = arg.info.fϕ
	ϕₘ   = arg.info.ϕₘ
	fstr = "  "*num2si(String,fₘ,1)*"Hz"
	ϕstr = num2si(String,ϕₘ,1)*"°  "

	link       := :x
	framestyle := [:semi :semi]
	layout     := (2,1)
	xscale     := :log10
	xlims      := (fₘᵢₙ,fₘₐₓ)
	seriestype := :line
	linewidth  --> 2
	legend     --> false
	grid       --> true
	
	@series begin
		subplot            := 1
		seriestype         := :line
		markershape        := :circle
		series_annotations := [(fstr, 10, :left, :bottom, 0.0)]
		[fₘ],[0]
	end
	@series begin
		subplot            := 2
		seriestype         := :line
		markershape        := :circle
		series_annotations := [(ϕstr, 10,  :right, :vcenter, 90.0)]
		[fₘ],[ϕₘ-180]
	end
	@series begin
		subplot := 1
		ylabel  := "Magnitude (dB)"
		xminorticks := true
		f, m
	end
	@series begin
		subplot := 2
		ylabel  := "Phase (°)"
		xlabel  := "Frequency (Hz)"
		xminorticks := true
		f, ϕ
	end
end




"""
	bodeinfo = pllbodeinfo(pll)

Return the phase margin and gain margin of a PLL as a [`BodeInfo`](@ref) object.
"""
function pllbodeinfo(pll::PLL)
	# To find the cutoff frequency and phase/gain margins, we need the bounds of the 
	# transfer function. The following is taken from _bounds_and_features in the 
	# ControlSystems package.
	OL = pll["Open-loop"]
	z = roots(OL.num)
	p = roots(OL.den)
	Tzp = promote_type(eltype(eltype(z)), eltype(eltype(p)))
	zp  = vcat(Tzp[], z..., p...)
	zp  = zp[imag(zp) .>= 0.0]
	ωlog = log10.(abs.(zp))
	ωflt = ωlog[ωlog .> -4]
	if !isempty(ωflt)
		favg = exp10(sum(ωflt)/length(ωflt))/(2*pi)
	else
		favg = 1e3
	end
	
	# Function for magnitude:
	m₀(f::Float64) = dB20(abs(OL(2im*pi*log10(f))))

	# Function for phase+π:
	ϕ₀(f::Float64) = angle(OL(2im*pi*log10(f)))+π

	if sign(m₀(favg)) != sign(m₀(eps(Float64)))
		f1 = eps(Float64)
		f2 = favg
	else
		f1 = favg
		f2 = 2*favg
		maxcount = 100
		while maxcount>0 && sign(m₀(f1)) != sign(m₀(f2))
			f2 *= 2
		end
		if maxcount==0
			fg0log = NaN
		else
			fg0log = find_zero(m₀, (f1,f2))
		end
	end
	# Find frequency of zero crossing:
	# fg0log = find_zero(m₀, favg, Order1())
	fϕ0log = find_zero(ϕ₀, favg, Order1())
	fg0 = 10^fg0log
	fϕ0 = 10^fϕ0log


	pm = 180*ϕ₀(fg0log)/pi
	gm = m₀(fϕ0log)

	PllBodeInfo(fg0,pm, fϕ0,gm)
end





"""
	pllbode(pll::PLL,f)

Calculate the frequency response of a PLL at frequencies in f.
"""
function pllbode end

pllbode(pll::PLL, f::Real) = pllbode(pll, [Float64(f)])
pllbode(pll::PLL, f::AbstractArray{<:Real}) = pllbode(pll, Float64.(f[:]))
function pllbode(pll::PLL, f::Vector{Float64})
	OL = pll["Open-loop"]
	m,ϕ = bode(OL,f)
	PllBodeResponse(f,m,ϕ,pllbodeinfo(pll))
end


function _bode_fvec(pll::PLL)
	fϕ = pllbodeinfo(pll).fϕ
	f1 = floor(log10(fϕ/1e3))
	f2 = ceil(log10(fϕ*1e3))
	10 .^LinRange(f1,f2,400)
end