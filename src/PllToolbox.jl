module PllToolbox

using Tables
using LinearAlgebra
using OrdinaryDiffEq
using Polynomials
using Printf
using Roots
using RecipesBase

"""
TODO!
1. Change how blocks are stored in the PLL. They should be accessible through their name, 
but also have a consistent ordering for proper NTF calculation.
2. Add function to return NTFs of all blocks, as a Dict probably. 
3. Add function that finds info on transfer functions, such as cross points, bandwidth,
type, etc.
4. Make analysis results accessible as DataFrames.
5. Add plot recipees for bode, noise, step.
6. Add function to calculate noise contribution by NTF type.
"""

include("types/tf.jl")
include("types/noise.jl")
include("types/block.jl")
include("types/other.jl")
include("types/pll.jl")
include("utilities/block_factory.jl")
include("utilities/filter_factory.jl")
include("utilities/num2si.jl")
include("analysis/step_response.jl")
include("analysis/frequency_response.jl")
include("analysis/noise.jl")




#-------------------------------------------------------------------------------------------
# PLL noise operations

export pllnoise

"""
	f,totnoise,blocknoise = pllnoise(pll,f)

Return the frequency vector, total noise and individual block noise of a PLL. 
Calculates the 
"""
# pllnoise(pll::PLL, f::Real) = pllnoise(pll, [Float64(f)])
# pllnoise(pll::PLL, f::AbstractArray{<:Real}) = pllnoise(pll, Float64.(f[:]))
# function pllnoise(pll::PLL, f::Vector{Float64})
# 	G,L,F,_,O = plltf(pll)
# 	nblkIp = _pllnoise_branch(pll.input, f, G*O)
# 	nblkFw = _pllnoise_branch(pll.forward, f, O/(1+L*F))
# 	nblkFb = _pllnoise_branch(pll.feedback, f, -O*G)
# 	nblkOp = _pllnoise_branch(pll.output, f)
# 	nblocks = [nblkIp;nblkFw;nblkFb;nblkOp]
# 	totnoise = zeros(size(f))
# 	for nb in nblocks 
# 		if ~isempty(nb.indnoise)
# 			totnoise = totnoise .+ sum(nb.indnoise,dims=2)
# 		end
# 	end
# 	(f, totnoise, nblocks)
# end




# function _pllnoise_branch(blocks::Vector{Block}, f::Vector{Float64}, h::LTISystem=tf(1))
# 	num_noisy = count(b->!isempty(b.noise), blocks)
# 	bnoise = Vector{BlockNoise}(undef, num_noisy)
# 	htot = tf(h)
# 	nᵢ = length(bnoise)
# 	for bᵢ in reverse(eachindex(blocks))
# 		if !isempty(blocks[bᵢ].noise)
# 			bnoise[nᵢ] = BlockNoise(blocks[bᵢ], f, htot)
# 			nᵢ -= 1
# 		end
# 		htot = htot*blocks[bᵢ].H
# 	end
# 	bnoise
# end




#-------------------------------------------------------------------------------------------
# PLL Bode operations






#-------------------------------------------------------------------------------------------
# Noise plot functions

export pllnoiseplot

"""
	pllnoiseplot(pll,f,interactive=false)
	pllnoiseplot(f,totalnoise,blocknoise,interactive=false)
"""
# pllnoiseplot(pll::PLL, f::AbstractArray{Float64}, interactive::Bool=false) = pllnoiseplot(pllnoise(pll,f)..., interactive)
# function pllnoiseplot(f::AbstractArray{Float64}, totalnoise::AbstractArray{Float64}, blocknoise::AbstractArray{BlockNoise}, interactive::Bool=false)
# 	indnoise = reduce(hcat, (n->n.indnoise).(blocknoise))
# 	srcnoise = reduce(hcat, (n->n.srcnoise).(blocknoise))
# 	noisestr = reduce(vcat, (n->n.noisestr).(blocknoise))
	
# 	if interactive
# 		_pllnoiseplot_interactive(f, totalnoise, indnoise, srcnoise, noisestr)
# 	else
# 		_pllnoiseplot_plain(f, totalnoise, indnoise, srcnoise, noisestr)
# 	end
# end


# function _pllnoiseplot_plain(f, totalnoise, indnoise, srcnoise, noisestr)
# 	yₘₐₓ = 10*log10(maximum(totalnoise))
# 	yₘᵢₙ = max(-200, 10*log10(minimum(hcat(indnoise,srcnoise))))

# 	Plots.plot(
# 		f,hcat(totalnoise,indnoise,srcnoise),
# 		label=vcat("Total",noisestr)
# 	)
# 	gui()
# end


# function _pllnoiseplot_interactive(f, totalnoise, indnoise, srcnoise, noisestr)

# 	yₘₐₓ = 10*log10(maximum(totalnoise))
# 	yₘᵢₙ = max(-200, 10*log10(minimum(hcat(indnoise,srcnoise))))

# 	line2d = matplotlib.lines.Line2D
# 	fig,ax = subplots(2,1,gridspec_kw=Dict("height_ratios"=>[1,40]))
# 	fig.set_visible(false)
# 	fig.canvas.manager.set_window_title("PLL phase noise")

# 	ltot = ax[2].semilogx(f, 10*log10.(totalnoise), lw=2, label="Total", c="black")
# 	lind = ax[2].semilogx(f, 10*log10.(indnoise), "--", label=noisestr)
# 	ax[2].set_prop_cycle(nothing)
# 	lsrc = ax[2].semilogx(f, 10*log10.(srcnoise), ":", label=noisestr, visible=false)

# 	ax[1].spines["left"].set_visible(false)
# 	ax[1].spines["right"].set_visible(false)
# 	ax[1].spines["top"].set_visible(false)
# 	ax[1].spines["bottom"].set_visible(false)
# 	ax[1].set_xticks([])
# 	ax[1].set_yticks([])
# 	ax[2].set_xlabel("Offset frequency (Hz)")
# 	ax[2].set_ylabel("PSD (dBc/Hz)")
# 	ax[2].set_xlim([minimum(f),maximum(f)])
# 	ax[2].set_ylim([yₘᵢₙ,yₘₐₓ+10])
# 	ax[2].grid(which="major", axis="both")
# 	ax[2].grid(which="minor", axis="both", linestyle=":")
# 	ax[2].xaxis.set_major_formatter(matplotlib.ticker.EngFormatter())
	
# 	line_legend = ax[2].legend(lind, noisestr)
# 	line_legend_texts = line_legend.get_texts()
# 	group_legend = ax[1].legend(
# 			fill(line2d([0],[0],ls="none"),4),
# 			["Total","Individual","Source","Legend"],
# 			loc="lower center", ncol=4, frameon=false, handletextpad=0
# 		)
# 	group_legend_texts = group_legend.get_texts()

# 	filter_type = Dict{eltype(line_legend_texts),String}()
# 	for txt in line_legend_texts
# 		filter_type[txt] = "line"
# 	end
# 	for txt in group_legend_texts
# 		filter_type[txt] = "group"
# 	end

# 	line_state  = Dict(txt=>true for txt in line_legend_texts)
# 	line_data   = Dict(zip(line_legend_texts,eachrow(hcat(lind,lsrc))))
	
# 	group_state = Dict(zip(group_legend_texts,[true,true,false,true]))
# 	group_check = Dict(zip(group_legend_texts,group_legend.get_lines()))
# 	group_data  = Dict(zip(group_legend_texts,[ltot,lind,lsrc,[line_legend]]))
	
# 	for txt in line_legend_texts
# 		txt.set_picker(true)
# 		if line_state[txt]
# 			txt.set_color("black")
# 		else
# 			txt.set_color("gray")
# 		end
# 	end

# 	for txt in group_legend_texts
# 		txt.set_picker(true)
# 		temp_p = group_check[txt].get_data()
# 		temp_x = sum(temp_p[1])/2
# 		temp_y = temp_p[2][1]
# 		group_check[txt].set_data(([temp_x],[temp_y]))
# 		group_check[txt].set_linestyle("none")
# 		group_check[txt].set_marker("o")
# 		group_check[txt].set_markeredgecolor("black")
# 		if group_state[txt]
# 			group_check[txt].set_markerfacecolor("black")
# 			txt.set_color("black")
# 		else
# 			group_check[txt].set_markerfacecolor("none")
# 			txt.set_color("gray")
# 		end
# 	end

# 	function _pllnoiseplot_callback(event)
# 		txt = event.artist
# 		if filter_type[txt]=="group"
# 			group_state[txt] = !group_state[txt]
# 			for d in group_data[txt]
# 				d.set_visible(group_state[txt])
# 			end
# 			if group_state[txt]
# 				group_check[txt].set_markerfacecolor("black")
# 				txt.set_color("black")
# 				for t in line_legend.get_texts()
# 					t.set_color("black")
# 				end
# 			else
# 				group_check[txt].set_markerfacecolor("none")
# 				txt.set_color("gray")
# 			end

# 		elseif filter_type[txt]=="line"
# 			line_state[txt] = !line_state[txt]
# 			line_data[txt][1].set_visible(line_state[txt]&&group_state[group_legend_texts[2]])
# 			line_data[txt][2].set_visible(line_state[txt]&&group_state[group_legend_texts[3]])
			
# 			if line_state[txt]
# 				txt.set_color("black")
# 			else
# 				txt.set_color("gray")
# 			end
# 		end
# 		draw()
# 	end

# 	fig.canvas.mpl_connect("pick_event", _pllnoiseplot_callback)
# 	fig.set_visible(true)
# 	display(fig)
# 	return fig
# end




export PllNoiseInfo
struct PllNoiseInfo
	fc::Float64
	fstart::Float64
	fstop::Float64
	ipn::Float64
	pe::Float64
	rfm::Float64
	j::Float64
end

function Base.show(io::IO, x::PllNoiseInfo)
	println("PLL phase noise info:")
	println("    center frequency: $(num2si(String,x.fc))Hz")
	println("    integration limits: $(num2si(String,x.fstart))Hz to $(num2si(String,x.fstop))Hz")
	println("    integrated phase noise: $(round(x.ipn,digits=2)) dBc")
	println("    phase error: $(round(x.pe,digits=2))°")
	println("    residual FM: $(num2si(String,x.rfm))Hz")
	println("    jitter: $(num2si(String, x.j))s")
end


export pllnoiseinfo
function pllnoiseinfo(pll::PLL, f::AbstractArray{Float64}, fc::Float64, fstart::Float64=f[1], fstop::Float64=f[end])
	f,pn, = pllnoise(pll,f)
	pllnoiseinfo(f,pn,fc,fstart,fstop)
end


function pllnoiseinfo(f::Vector{Float64}, pn::AbstractArray{Float64}, fc::Float64, fstart::Float64=f[1], fstop::Float64=f[end])
	istart = findfirst(x->x>=fstart, f)
	istop = findlast(x->x<=fstop, f)
	fi  = view(f, istart:istop)
	pni = view(pn, istart:istop,:)
	fd  = view(f, istart:istop-1)
	pnd = view(pn, istart:istop-1,:)
	df = diff(fi)

	ipn = 2*sum(pnd .* df)
	pe  = 180*sqrt(ipn)/π
	rfm = sqrt(2*sum(fd.^2 .* pnd .* df))
	j   = pe./(360*fc)

	return PllNoiseInfo(fc,fstart,fstop,ipn,pe,rfm,j)
end


export pllnoiseprint
function pllnoiseprint(pll::PLL, f::AbstractArray{Float64}, fc::Float64, fstart::Float64=f[1], fstop::Float64=f[end])
	show(pllnoiseinfo(pll,f,fc,fstart,fstop))
end

end # module
