module PllToolbox

using Tables
using LinearAlgebra
using OrdinaryDiffEq
using Polynomials
using Printf
using Roots
using RecipesBase


export pllnoise, pllnoiseinfo, pllintegnoise, noiseplot, noiseplot!, contributionplot, contributionplot!
export pllstep, stepplot, stepplot!
export dB10, dB20, PLL, pllntf
export pllbode, pllbodeinfo
export Block
export TF, bode, unwrap, unwrap!, classify
export AbstractNoise, NoNoise, WhiteNoise, PinkNoise, ΣΔNoise, ΔΣNoise, SDNoise, DSNoise
export pllref, pllgain, plldiv, pllmmd, pllvco
export pllfilter

include("types/tf.jl")
include("types/noise.jl")
include("types/block.jl")
include("types/pll.jl")
include("utilities/block_factory.jl")
include("utilities/filter_factory.jl")
include("utilities/num2si.jl")
include("analysis/step_response.jl")
include("analysis/frequency_response.jl")
include("analysis/noise.jl")
include("documentation.jl")
end # module
