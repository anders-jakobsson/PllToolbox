
#-------------------------------------------------------------------------------------------
# Reference clock source

export pllref
"""
	pllref(name::String, noise::Vector{AbstractNoise}=[], descr::String="")

Create a reference clock with the given name and noise source. `noise` can be either an 
instance of `Noise` or a vector of `Noise` instances.
"""
function pllref end
pllref(name::String, noise::Vector{<:AbstractNoise}=AbstractNoise[], descr::String="") = Block(name, TF(1,1,name), noise, descr)
pllref(name::String, noise::AbstractNoise, descr::String="") = pllref(name, [noise], descr)




#-------------------------------------------------------------------------------------------
# Gain block

export pllgain
"""
	pllgain(name::String, a::Number=1, noise::Vector{AbstractNosie}=[], descr::String="")

Create a gain block with the given name, gain and noise source. `noise` can be either an 
instance of [`Noise`](@ref) or a vector of [`Noise`](@ref) instances.
"""
function pllgain end
pllgain(name::String, a::Number=1, noise::Vector{<:AbstractNoise}=AbstractNoise[], descr::String="") = Block(name, TF(a,1,name), noise, descr)
pllgain(name::String, a::Number, noise::AbstractNoise, descr::String="") = pllgain(name, a, [noise], descr)




#-------------------------------------------------------------------------------------------
# Frequency divider block

export plldiv
"""
	plldiv(name::String, ndiv::Int=2, noise::Vector{AbstractNoise}=[], descr::String="")

Create an integer divider with the given name, divider ratio and noise source. `noise` can 
be either an instance of `Noise` or a vector of `Noise` instances.
"""
function plldiv end
plldiv(name::String, ndiv::Int=2, noise::Vector{<:AbstractNoise}=AbstractNoise[], descr::String="") = Block(name, TF(1,ndiv,name), noise, descr)
plldiv(name::String, ndiv::Int, noise::AbstractNoise, descr::String="") = plldiv(name, ndiv, [noise], descr)




#-------------------------------------------------------------------------------------------
# Multi-modulus divider (MMD) block

export pllmmd
"""
	pllmmd(name::String, ndiv::Real=100, fs::Real, order::Int=3, noise::Vector{AbstractNoise}=[], descr::String="")

Create a multi-modulus divider (MMD) with ΣΔ-modulation. A ΣΔ-noise instance is 
automatically created and added to the block. `noise` can be either an instance of 
[`Noise`](@ref) or a vector of [`Noise`](@ref) instances.
"""
function pllmmd end

function pllmmd(name::String, ndiv::Real, Fs::Real, order::Int=3, noise::Vector{<:AbstractNoise}=AbstractNoise[], descr::String="")
	ΣΔn = ΣΔNoise(name*" ΣΔ-noise", (π/ndiv)^2/(3Fs), 2, 2order-2, Fs/pi)
	Block(name, TF(1,ndiv,name), [ΣΔn;noise], descr)
end
pllmmd(name::String, ndiv::Real, Fs::Real, order::Int, noise::AbstractNoise, descr::String="") = pllmmd(name, ndiv, Fs, order, [noise], descr)




#-------------------------------------------------------------------------------------------
# Voltage controlled oscillator (VCO) block

export pllvco
"""
	pllvco(name::String, kvco::Real, noise::Vector{<:AbstractNoise}=AbstractNoise[], descr::String="")

Create a VCO block with the given name, gain in Hz/V and noise. 
"""
function pllvco end

pllvco(name::String, kvco::Real, noise::Vector{<:AbstractNoise}=AbstractNoise[], descr::String="") = Block(name, TF(kvco*2π,[1,0],name), noise, descr)
pllvco(name::String, kvco::Real, noise::AbstractNoise, descr::String="") = pllvco(name, kvco, [noise], descr)
