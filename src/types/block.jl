using AbstractTrees




#------------------------------------------------------------------------------------------
# Block type
export Block
"""
    Block

Representation of a block in the PLL model. Blocks are normally constructed using one of the 
factory functions, see [`pllref`](@ref), [`pllgain`](@ref), [`plldiv`](@ref),
[`pllmmd`](@ref) and [`pllvco`](@ref).
"""
struct Block{T}
	name::String
	tf::TF
	noise::Vector{T}
	descr::String

	@doc """
		Block(name, tf::TF=TF(1), noise::Vector{<:AbstractNoise}, description::String="")

	Construct a PLL block from the given name, transfer function and noise argument.
	The name must be a non-empty string. See also [`pllref`](@ref), [`pllgain`](@ref), 
	[`plldiv`](@ref), [`pllmmd`](@ref) and [`pllvco`](@ref).


	# Examples
	```@meta
	DocTestSetup = quote
		using PllToolbox
	end
	```
	```jldoctest
	julia> Block("VCO", TF(1,[1,0]), nothing)

	julia> Block("")

	"""
	function Block(name::String, tf::TF=TF(name,1), noise::Vector{T}=AbstractNoise[], descr::String="") where {T<:AbstractNoise}
		# Hconv = convert(TransferFunction{Continuous,ControlSystems.SisoRational{Float64}}, H)
		new{T}(name, tf, noise, descr)
	end
end

Block(name::String, tf::TF, noise::AbstractNoise, descr::String="") = Block(name, tf, [noise], descr)




"""
    Block(b, name=b.name, tf=b.tf, noise=b.noise)

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
Block(block::Block; name::String=block.name, tf::TF=block.tf, noise::AbstractNoise=block.noise, descr::String=block.descr) = Block(name, tf, noise, descr)




# Custom pretty-printing.
function Base.show(io::IO, x::Block)
	name = x.name
	if isnothing(x.noise)
		print(io, "PLL block \"$(x.name)\" with no noise sources")
	elseif length(x.noise)>1
		print(io, "PLL block \"$(x.name)\" with $(length(x.noise)) noise sources")
	else
		print(io, "PLL block \"$(x.name)\" with one noise source")
	end
end








#-------------------------------------------------------------------------------------------
# CascadedBlock type
struct CascadedBlock
	block::Block
	ntf::TF
end


Base.propertynames(cb::CascadedBlock, private::Bool=false) = (fieldnames(Block)...,:ntf)
function Base.getproperty(cb::CascadedBlock, name::Symbol)
	return name===:ntf ? getfield(cb,:ntf) : getfield(getfield(cb,:block),name)
end


# Custom pretty-printing.
Base.show(io::IO, x::CascadedBlock) = show(io, getfield(x,:block))






AbstractTrees.children(node::CascadedBlock) = ()