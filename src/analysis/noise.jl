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


struct PllIntegNoise{T}
	frequency::Pair{T,T}
	total::T
	byname::Dict{String,T}
end


struct PllNoiseInfo
	fc::Float64
	fstart::Float64
	fstop::Float64
	ipn::Float64
	pe::Float64
	rfm::Float64
	j::Float64
end


pllnoise(pll::PLL, f::Real, kw...) = pllnoise(pll, [Float64(f)]; kw...)
pllnoise(pll::PLL, f::AbstractArray{<:Real}, kw...) = pllnoise(pll, Float64.(f[:]); kw...)
function pllnoise(pll::PLL, f::Vector{Float64}; byfiltertype::Bool=false)
	allblocks = vcat(pll.input,pll.forward,pll.feedback,pll.output)
	if byfiltertype
		numnoise = length(unique(block.ntf.type for block in allblocks))
	else
		numnoise = sum([length(b.noise) for b in allblocks])
	end
	
	total  = zeros(length(f))
	byname = Dict{String,PllNoiseSource{Float64}}()
	for bᵢ in allblocks 
		pr = magresp(bᵢ.ntf,f).^2
		for nᵢ in bᵢ.noise
			noise = nᵢ(f)
			contr = noise.*pr
			if byfiltertype
				type = bᵢ.ntf.type
				if String(type) in keys(byname)
					nsrc = byname[String(type)]
					nsrc.source .+= noise
					nsrc.contribution .+= contr
				else
					byname[String(type)] = PllNoiseSource{Float64}(type, noise, contr)
				end
			else
				byname[nᵢ.name] = PllNoiseSource{Float64}(bᵢ.ntf.type,noise,contr)
			end
			total .+= contr
		end
	end

	PllNoise(f, total, collect(values(byname)), byname)
end


function pllintegnoise(pll::PLL, fstart::Float64=1e3, fstop::Float64=100e6, npoints::Integer=1000; kw...)
	f = 10 .^LinRange(log10(fstart),log10(fstop),npoints)
	pllintegnoise(pll, f; kw...)
end
function pllintegnoise(pll::PLL, f::Vector{Float64}; kw...)
	pllintegnoise(pllnoise(pll,f;kw...))
end
function pllintegnoise(pn::PllNoise, fstart::Number=pn.frequency[1], fstop::Number=pn.frequency[end])
	istart = findfirst(x->x>=fstart, pn.frequency)
	istop  = findlast(x->x<=fstop, pn.frequency)
	df = diff(pn.frequency[istart:istop])

	total  = 2*sum(pn.total[istart:istop-1] .* df)
	byname = Dict{String,Float64}()
	for (name,pns) in pn.byname
		contr = 2*sum(pns.contribution[istart:istop-1] .* df)
		byname[name] = contr
	end
	PllIntegNoise(Pair(pn.frequency[istart],pn.frequency[istop]), total, byname)
end





function pllnoiseinfo(pll::PLL, fc::Float64, fstart::Float64=1e3, fstop::Float64=100e6, npoints::Integer=1000)
	pn = pllnoise(pll,10 .^LinRange(log10(fstart),log10(fstop),npoints))
	pllnoiseinfo(pn,fc)
end


function pllnoiseinfo(pn::PllNoise, fc::Float64, fstart::Float64=pn.frequency[1], fstop::Float64=pn.frequency[end])
	f = pn.frequency
	istart = findfirst(x->x>=fstart, f)
	istop = findlast(x->x<=fstop, f)
	fi  = view(f, istart:istop)
	fd  = view(f, istart:istop-1)
	pnd = view(pn.total, istart:istop-1,:)
	df = diff(fi)

	ipn = 2*sum(pnd .* df)
	pe  = 180*sqrt(ipn)/π
	rfm = sqrt(2*sum(fd.^2 .* pnd .* df))
	j   = pe./(360*fc)

	return PllNoiseInfo(fc,fstart,fstop,ipn,pe,rfm,j)
end




#-------------------------------------------------------------------------------------------
# Show implementation
function Base.show(io::IO, ::MIME"text/plain", x::PllNoiseInfo)
	println("PLL phase noise info:")
	println("    center frequency: $(num2si(String,x.fc))Hz")
	println("    integration limits: $(num2si(String,x.fstart))Hz to $(num2si(String,x.fstop))Hz")
	println("    integrated phase noise: $(num2si(String,x.ipn)) dBc")
	println("    phase error: $(num2si(String,x.pe))°")
	println("    residual FM: $(num2si(String,x.rfm))Hz")
	println("    jitter: $(num2si(String, x.j))s")
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
	title       --> "PLL phase noise"
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




@userplot IntegNoisePlot

@recipe function finp(cp::IntegNoisePlot)
	length(cp.args)≥1 || error("No arguments supplied")
	arg  = cp.args[1]
	f    = arg.frequency
	tot  = arg.total
	title  --> "Noise contribution"
	@series begin
		seriestype := :pie
		x = String.(keys(arg.byname))
		y = collect(values(arg.byname))
		x,y	
	end
end



@userplot ContributionPlot

contributionplot(pll::PLL; kw...) = contributionplot(pll,_noise_fvec(pll); kw...)
contributionplot(pll::PLL,f::AbstractVector{T}; kw...) where T<:Number = contributionplot(pllnoise(pll,f); kw...)
contributionplot!(pll::PLL; kw...) = contributionplot!(pll,_noise_fvec(pll); kw...)
contributionplot!(pll::PLL,f::AbstractVector{T}; kw...) where T<:Number = contributionplot!(pllnoise(pll,f); kw...)

@recipe function fcp(cp::ContributionPlot)
	length(cp.args)≥1 || error("No arguments supplied")
	arg  = cp.args[1]
	f    = arg.frequency
	tot  = arg.total
	fₘᵢₙ = minimum(f)
	fₘₐₓ = maximum(f)
	
	title --> "Noise contribution"
	if length(f)>1
		framestyle  := :semi
		xlims       := (fₘᵢₙ,fₘₐₓ)
		ylims       := (0,100)
		xminorticks := true
		xscale      := :log10
		xlabel      := "Offset frequency (Hz)"
		ylabel      := "Noise contribution (%)"
		seriestype  := :line
		markershape := :none
		legend     --> true
		grid       --> false
		
		fill = zeros(length(f))
		for (name,src) in arg.byname
			contr = 100*src.contribution./tot
			y = fill .+ contr
			@series begin
				label     := name
				fillrange := fill
				f, y
			end
			fill = y
		end
	else
		@series begin
			seriestype := :pie
			x = String.(keys(arg.byname))
			y = [src.contribution[1] for src in values(arg.byname)]
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