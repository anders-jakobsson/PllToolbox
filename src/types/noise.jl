
# ------------------------------------------------------------------------------------------
# Imported types

using Statistics:std
using Interpolations:BSpline,Gridded,Line,Linear,OnGrid,Quadratic,extrapolate,interpolate,scale



#-------------------------------------------------------------------------------------------
# Declared types

export AbstractNoise, NoNoise, WhiteNoise, PinkNoise
export ΣΔNoise, ΔΣNoise, SDNoise, DSNoise




"""
    AbstractNoise

Abstract functor supertype for representing a PLL noise source.

A sub-type of AbstractNoise must implement the following fields and methods:

    name::String                          name of noise source
    (<:AbstractNoise)(Vector{Float64})    functor method that calculates PSD
"""
abstract type AbstractNoise end




"""
    (n::<:AbstractNoise)(f::Vector{Float64})

Return the noise power spectral density (PSD) at each frequency in f.
"""
function (n::AbstractNoise)(f::Vector{Float64}) end








"""
	WhiteNoise(name::String, psd::Float64, H::LTISystem=tf(1))

Construct a white noise source from the given name, power spectral density and transfer 
function. `H` is an additional transfer function applied when calculating the PSD.

# Examples
```@meta
DocTestSetup = quote
    using PllToolbox
end
```
```jldoctest
julia> WhiteNoise("R1 noise", 4.0*1.380649e-23*300*1e3) # Thermal noise of 1kohm resistor

"""
struct WhiteNoise <: AbstractNoise
	name::String
	psd::Float64
	tf::TF
	function WhiteNoise(name::String, psd::Real, tf::TF=TF())
		new(name, psd, tf)
	end
end


function (x::WhiteNoise)(f::Vector{Float64})
	m = abs.(magresp(x.tf,f))
	x.psd * (m.^2)
end








"""
	PinkNoise(name, fx, px, logy=true, H=tf(1))

Construct a pink noise source, interpolated from a set of frequency/PSD pairs.

# Arguments
- `name::String`: noise source name
- `fx::Vector{<:Real}`: vector of frequencies
- `px::Vector{<:Real}`: vector of power spectral densities
- `logy::Bool`: set to `true` if px is in dB/Hz
- `H::LTISystem`: additional transfer function to apply when calculating the PSD


# Examples
```@meta
DocTestSetup = quote
    using PllToolbox
endf
```
```jldoctest
julia> PinkNoise("Colored noise", 10 .^(2:6), -[100,120,135,150,160])

"""
struct PinkNoise <: AbstractNoise
	name::String
	fx::Vector{Float64}
	px::Vector{Float64}
	logy::Bool
	tf::TF
	function PinkNoise(name::AbstractString, fx::AbstractVector{<:Real}, px::AbstractVector{<:Real}, logy::Bool=true, tf::TF=TF())
		if length(fx)<2
			error("frequency vector must contain at least two elements")
		end
		if length(fx)!=length(px)
			error("frequency and PSD vectors must have the same length")
		end

		if !logy
			px = 10*log10.(px[:])
		end
		new(name, fx[:], px[:], logy, tf)
	end
end

function (x::PinkNoise)(f::Vector{Float64})
	m = abs.(magresp(x.tf,f))
	y = similar(f)

	dflin = diff(x.fx)
	dflog = diff(log10.(x.fx))
	if isnothing(findfirst(x->abs(x-dflog[1])>eps(0.0), dflog))
		frange = range(log10(x.fx[1]), log10(x.fx[end]), length(x.fx))
		itp = scale(interpolate(x.px, BSpline(Quadratic(Line(OnGrid())))), frange)
		etp = extrapolate(itp,Line())
		y = etp(log10.(f))

	elseif isnothing(findfirst(x->abs(x-dflin[1])>eps(0.0), dflin))
		frange = range(x.fx[1], x.fx[end], length(x.fx))
		itp = scale(interpolate(x.px, BSpline(Quadratic(Line(OnGrid())))), frange)
		etp = extrapolate(itp,Line())
		y = etp(f)
	else
		if std(dflin)<std(dflog)
			itp = interpolate((x.fx,), x.px, Gridded(Linear()))
			etp = extrapolate(itp,Line())
			y = etp(f)
		else
			itp = interpolate((log10.(x.fx),), x.px, Gridded(Linear()))
			etp = extrapolate(itp,Line())
			y = etp(log10.(f))
		end
	end

	return 10 .^(y[:]/10) .* (m.^2)
end








"""
	ΣΔNoise(name::String, a::Real, b::Real, c::Real, f0::Real, H::LTISystem=tf(1))
	ΔΣNoise(...)
	SDNoise(...)
	DSNoise(...)

Construct a ΣΔ-modulated noise source with the PSD a⋅[b·sin(f/f0)]^c. This type has several
aliases for convenience. `H` is an additional transfer function applied when calculating 
the PSD.
 """
struct ΣΔNoise <: AbstractNoise
	name::String
	a::Float64
	b::Float64
	c::Float64
	f0::Float64
	tf::TF
	function ΣΔNoise(name::String, a::Real, b::Real, c::Real, f0::Real, tf::TF=TF())
		new(name, a, b, c, f0, tf)
	end
end

function (x::ΣΔNoise)(f::Vector{Float64})
	m = abs.(magresp(x.tf,f))
	(x.a*(x.b*sin.(f/x.f0)).^x.c) .* (m.^2)
end


ΔΣNoise = ΣΔNoise
SDNoise = ΣΔNoise
DSNoise = ΣΔNoise
