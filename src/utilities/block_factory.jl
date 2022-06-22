
#--------------------------------------------------------------------------------------------------
# Reference clock source

export pllref
pllref(name::String, noise::Vector{<:AbstractNoise}=AbstractNoise[]) = Block(name, tf(1), noise)
pllref(name::String, noise::AbstractNoise) = pllref(name, [noise])




#--------------------------------------------------------------------------------------------------
# Frequency divider block

export plldiv
plldiv(name::String, ndiv::Int=2, noise::Vector{<:AbstractNoise}=AbstractNoise[]) = Block(name, tf(1,[ndiv]), noise)
plldiv(name::String, ndiv::Int, noise::AbstractNoise) = plldiv(name, ndiv, [noise])




#--------------------------------------------------------------------------------------------------
# Gain block

export pllgain
pllgain(name::String, a::Number=1, noise::Vector{<:AbstractNoise}=AbstractNoise[]) = Block(name, tf(a), noise)
pllgain(name::String, a::Number, noise::AbstractNoise) = pllgain(name, a, [noise])




#--------------------------------------------------------------------------------------------------
# Multi-modulus divider (MMD) block

export pllmmd
function pllmmd(name::String, Fs::Real, ndiv::Real=100, order::Int=3, noise::Vector{<:AbstractNoise}=AbstractNoise[])
	ΣΔn = ΣΔNoise(name*" quantization noise", (π/ndiv)^2/(3Fs), 2, 2order-2, Fs/pi)
	Block(name, tf(1,[ndiv]), [ΣΔn;noise])
end
pllmmd(name::String, Fs::Real, ndiv::Real, order::Int, noise::AbstractNoise) = pllmmd(name, Fs, ndiv, order, [noise])




#--------------------------------------------------------------------------------------------------
# Voltage controlled oscillator (VCO) block

export pllvco
pllvco(name::String, kvco::Real, noise::Vector{<:AbstractNoise}=AbstractNoise[]) = Block(name, tf(kvco*2π,[1,0]), noise)
pllvco(name::String, kvco::Real, noise::AbstractNoise) = pllvco(name, kvco, [noise])
