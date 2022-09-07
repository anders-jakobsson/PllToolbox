
# -----------------------------------------------------------------------------
# Imported types

using ControlSystems:LTISystem




#------------------------------------------------------------------------------
# Declared types
struct BlockNoise
	block::Block
	ntf::LTISystem
	srcnoise::Array{Float64,2}
	indnoise::Array{Float64,2}
	noisestr::Vector{String}
	function BlockNoise(block::Block, f::Vector{Float64}, ntf::LTISystem)
		magresp = abs.(freqresp(ntf, 2π*f)[:]).^2
		srcnoise = Array{Float64,2}(undef,length(f),length(block.noise))
		indnoise = Array{Float64,2}(undef,length(f),length(block.noise))
		noisestr = Vector{String}(undef,length(block.noise))
		for i in eachindex(block.noise)
			srcnoise[:,i] = block.noise[i](f)
			indnoise[:,i] = srcnoise[:,i].*magresp
			noisestr[i]   = block.noise[i].name
		end
		new(block,ntf,srcnoise,indnoise,noisestr)
	end
end




#-------------------------------------------------------------------------------------------
# Interface implementations
function Base.show(io::IO, x::BlockNoise)
	if length(x.noisestr)>1
		numstr = "$(length(x.noisestr)) noise sources"
	elseif length(x.noisestr)>0
		numstr = "one noise source"
	else
		numstr = "no noise"
	end
	if length(num(x.ntf)[1])>1 && length(den(x.ntf)[1])>1
		ntfstr = "NTF:\n" * plltfrepr(x.ntf)
	else
		ntfstr = "noise gain: $(num(x.ntf)[1][1]/den(x.ntf)[1][1])"
	end
	println(io, "PLL block $(x.block.name) with $(numstr), $(ntfstr)")
	return
end
