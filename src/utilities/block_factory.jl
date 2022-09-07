
#--------------------------------------------------------------------------------------------------
# Reference clock source

export pllref
"""
	pllref(name, noise=[])

Create a reference clock with the given name and noise source. `noise` can be either an 
instance of `Noise` or a vector of `Noise` instances.
"""
function pllref(name, noise) end
pllref(name::String, noise::Vector{<:AbstractNoise}=AbstractNoise[]) = Block(name, tf(1), noise)
pllref(name::String, noise::AbstractNoise) = pllref(name, [noise])




#--------------------------------------------------------------------------------------------------
# Frequency divider block

export plldiv
"""
	plldiv(name, ndiv=2, noise=[])

Create an integer divider with the given name, divider ratio and noise source. `noise` can 
be either an instance of `Noise` or a vector of `Noise` instances.
"""
function plldiv(name, ndiv, noise) end
plldiv(name::String, ndiv::Int=2, noise::Vector{<:AbstractNoise}=AbstractNoise[]) = Block(name, tf(1,[ndiv]), noise)
plldiv(name::String, ndiv::Int, noise::AbstractNoise) = plldiv(name, ndiv, [noise])




#--------------------------------------------------------------------------------------------------
# Gain block

export pllgain
"""
	pllgain(name, a=1, noise=[])

Create a gain block with the given name, gain and noise source. `noise` can be either an 
instance of [`Noise`](@ref) or a vector of [`Noise`](@ref) instances.
"""
function pllgain(name, a, noise) end
pllgain(name::String, a::Number=1, noise::Vector{<:AbstractNoise}=AbstractNoise[]) = Block(name, tf(a), noise)
pllgain(name::String, a::Number, noise::AbstractNoise) = pllgain(name, a, [noise])




#-------------------------------------------------------------------------------------------
# Multi-modulus divider (MMD) block

export pllmmd
"""
	pllmmd(name, ndiv=100, fs, order=3, noise=[])

Create a multi-modulus divider (MMD) with ΣΔ-modulation. A ΣΔ-noise instance is 
automatically created and added to the block. 
and noise source. `noise` can be either an instance of [`Noise`](@ref) or a vector of 
[`Noise`](@ref) instances.
"""
function pllmmd(name, Fs, ndiv, order, noise) end

function pllmmd(name::String, Fs::Real, ndiv::Real=100, order::Int=3, noise::Vector{<:AbstractNoise}=AbstractNoise[])
	ΣΔn = ΣΔNoise(name*" quantization noise", (π/ndiv)^2/(3Fs), 2, 2order-2, Fs/pi)
	Block(name, tf(1,[ndiv]), [ΣΔn;noise])
end
pllmmd(name::String, Fs::Real, ndiv::Real, order::Int, noise::AbstractNoise) = pllmmd(name, Fs, ndiv, order, [noise])




#-------------------------------------------------------------------------------------------
# Voltage controlled oscillator (VCO) block

export pllvco
"""
	pllvco(name::String, kvco::Real, noise::Vector{<:AbstractNoise}=AbstractNoise[])

Create a VCO block with the given name, gain in Hz/V and noise. 
"""
function pllvco(name, kvco, noise=[]) end

pllvco(name::String, kvco::Real, noise::Vector{<:AbstractNoise}=AbstractNoise[]) = Block(name, tf(kvco*2π,[1,0]), noise)
pllvco(name::String, kvco::Real, noise::AbstractNoise) = pllvco(name, kvco, [noise])
