# -----------------------------------------------------------------------------------------
# Imported types

using ControlSystems
using AbstractTrees




#------------------------------------------------------------------------------------------
# Declared types

export Block




"""
    Block(name, H=tf(1), noise=Vector{<:
    	AbstractNoise}())

Construct a PLL block from the given name, transfer function and noise argument.

The name must be a non-empty string. Blocks are normally constructed by one of the 
factory functions, such as pllgain or pllvco.


# Examples
```@meta
DocTestSetup = quote
    using PllToolbox
end
```
```jldoctest
julia> Block("VCO", tf(1,[1,0]), nothing)

julia> Block("")

"""
struct Block{T}
	name::String
	H::LTISystem
	noise::Vector{T}

	function Block(name::String, H::LTISystem=tf(1), noise::Vector{T}=AbstractNoise[]) where {T<:AbstractNoise}
		# Hconv = convert(TransferFunction{Continuous,ControlSystems.SisoRational{Float64}}, H)
		new{T}(name, H, noise)
	end
end

Block(name::String, H::LTISystem, noise::AbstractNoise) = Block(name, H, [noise])








"""
    Block(b, name=b.name, H=b.H, noise=b.noise)

Copy a block, with optional override of attributes.


# Examples
```@meta
DocTestSetup = quote
    using PllToolbox
end
```
```jldoctest
julia> b = Block("VCO", tf(1,[1,0]), nothing)

julia> Block(b, name="VCO2")

"""
Block(block::Block, name::String=block.name, H::LTISystem=block.tf, noise::AbstractNoise=block.noise) = Block(name, H, noise)








# Custom pretty-printing.
function Base.show(io::IO, x::Block)
	name = x.name
	if isnothing(x.noise)
		println(io, "PLL block \"$(x.name)\" with no noise sources")
	elseif length(x.noise)>1
		println(io, "PLL block \"$(x.name)\" with $(length(x.noise)) noise sources")
	else
		println(io, "PLL block \"$(x.name)\" with one noise source")
	end
end








AbstractTrees.children(node::Block) = ()