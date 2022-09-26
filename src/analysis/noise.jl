struct PllNoiseSource{T}
	type::Symbol
	source::Vector{T}
	contribution::Vector{T}
end


struct PllNoise{T}
	frequency::Vector{T}
	total::Vector{T}
	byindex::Vector{PllNoiseSource{T}}
	byname::Dict{String,PllNoiseSource{T}}
end


pllnoise(pll::PLL, f::Real) = pllnoise(pll, [Float64(f)])
pllnoise(pll::PLL, f::AbstractArray{<:Real}) = pllnoise(pll, Float64.(f[:]))
function pllnoise(pll::PLL, f::Vector{Float64})
	allblocks = vcat(pll.input,pll.forward,pll.feedback,pll.output)
	numnoise  = sum([length(b.noise) for b in allblocks])

	total   = zeros(length(f))
	byindex = Vector{PllNoiseSource{Float64}}(undef,numnoise)
	byname  = Dict{String,PllNoiseSource{Float64}}()
	kᵢ = 1
	for bᵢ in allblocks 
		pr = magresp(bᵢ.ntf,f).^2
		for nᵢ in bᵢ.noise
			noise = nᵢ(f)
			contr = noise.*pr
			nsrc  = PllNoiseSource{Float64}(bᵢ.ntf.type,noise,contr)
			byindex[kᵢ] = nsrc
			byname[nᵢ.name] = nsrc
			total .+= contr
			kᵢ += 1
		end
	end

	PllNoise(f, total, byindex, byname)
end




#-------------------------------------------------------------------------------------------
# Tables implementation

struct PllNoiseWrapper{T} <: Tables.AbstractColumns
	pn::PllNoise{T}
end

Tables.istable(::PllNoise) = true
Tables.columnaccess(::PllNoise) = true
Tables.columns(x::PllNoise) = PllNoiseWrapper(x)
function Tables.columnnames(x::PllNoiseWrapper)
	pn = getfield(x,:pn)
	Symbol.(vcat("Frequency","Total",collect(keys(pn.byname))))
end
function Tables.getcolumn(x::PllNoiseWrapper, k::Int)
	pn = getfield(x,:pn)
	k==1 && return pn.frequency
	k==2 && return pn.total
	pn.byindex[k-2].contribution
end
function Tables.getcolumn(x::PllNoiseWrapper, name::Symbol)
	pn = getfield(x,:pn)
	name==:Frequency && return pn.frequency
	name==:Total && return pn.total
	return pn.byname[String(name)].contribution
end
function Tables.schema(x::PllNoiseWrapper{T}) where T
	pn = getfield(x,:pn)
	Tables.Schema(Tables.columnnames(x),fill(T,2+length(pn.byindex)))
end




#-------------------------------------------------------------------------------------------
# Plot recipees

@userplot NoisePlot

noiseplot(pll::PLL; kw...) = noiseplot(pll,_noise_fvec(pll); kw...)
noiseplot(pll::PLL,f::AbstractVector{T}; kw...) where T<:Number = noiseplot(pllnoise(pll,f); kw...)
noiseplot!(pll::PLL; kw...) = noiseplot!(pll,_noise_fvec(pll); kw...)
noiseplot!(pll::PLL,f::AbstractVector{T}; kw...) where T<:Number = noiseplot!(pllnoise(pll,f); kw...)
	
@recipe function fnp(x::NoisePlot)
	length(x.args)≥1 || error("No arguments supplied")
	arg  = x.args[1]
	f    = arg.frequency
	tot  = dB10.(arg.total)
	fₘᵢₙ = minimum(f)
	fₘₐₓ = maximum(f)
	pₘᵢₙ = minimum(tot)
	pₘₐₓ = maximum(tot)
	
	width = get(plotattributes, :linewidth, 1)
	framestyle  := :semi
	xlims       := (fₘᵢₙ,fₘₐₓ)
	ylims       := (pₘᵢₙ-10,pₘₐₓ+10)
	xscale      := :log10
	seriestype  := :path
	xminorticks := true
	xlabel      := "Offset frequency (Hz)"
	ylabel      := "Phase noise (dBc/Hz)"
	legend      --> true
	grid        --> true
	
	@series begin
		linewidth := 2*width
		label     := "Total"
		f, tot
	end
	linewidth := width
	linestyle := :dash
			
	for (name,src) in arg.byname
		@series begin
			label := name
			f, dB10.(src.contribution)
		end
	end
end


@userplot ContributionPlot

contributionplot(pll::PLL; kw...) = contributionplot(pll,_noise_fvec(pll); kw...)
contributionplot(pll::PLL,f::AbstractVector{T}; kw...) where T<:Number = contributionplot(pllnoise(pll,f); kw...)
contributionplot!(pll::PLL; kw...) = contributionplot!(pll,_noise_fvec(pll); kw...)
contributionplot!(pll::PLL,f::AbstractVector{T}; kw...) where T<:Number = contributionplot!(pllnoise(pll,f); kw...)

@recipe function fcp(cp::ContributionPlot; byfiltertype=false)
	length(cp.args)≥1 || error("No arguments supplied")
	arg  = cp.args[1]
	f    = arg.frequency
	tot  = arg.total
	fₘᵢₙ = minimum(f)
	fₘₐₓ = maximum(f)
	
	group = Dict{String,Vector{Float64}}()
	if byfiltertype
		names = unique([String(src.type) for src in arg.byindex])
		for src in arg.byindex
			group[String(src.type)] = src.contribution .+ get(group,String(src.type),zeros(length(f)))
		end
	else
		for (name,src) in arg.byname
			group[name] = src.contribution
		end
		names = keys(group)
	end

	if length(f)>1
		framestyle  := :semi
		xlims       := (fₘᵢₙ,fₘₐₓ)
		ylims       := (0,100)
		xminorticks := true
		xscale      := :log10
		xlabel      := "Offset frequency (Hz)"
		ylabel      := "Noise contribution (%)"
		legend     --> true
		grid       --> false
		seriestype  := :line
		markershape := :none
		
		fill = zeros(length(f))
		for (name,n) in group
			contr = 100*n./tot
			y = fill .+ contr
			@series begin
				label       := name
				fillrange   := fill
				f, y
			end
			fill = y
		end
	else
		@series begin
			seriestype := :pie
			x = String.(keys(group))
			y = [g[1] for g in values(group)]
			x,y	
		end
	end
end


function _noise_fvec(pll::PLL)
	fϕ = pllbodeinfo(pll).fϕ
	f1 = floor(log10(fϕ/1e3))
	f2 = ceil(log10(fϕ*1e3))
	10 .^LinRange(f1,f2,400)
end