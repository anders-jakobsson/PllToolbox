export pllbodeplot
"""
	pllbodeplot(pll,f)

Create a Bode plot for a PLL at frequencies in f.
"""
pllbodeplot(pll::PLL, f::Real) = pllbodeplot(pll, [Float64(f)])
pllbodeplot(pll::PLL, f::AbstractArray{<:Real}) = pllbodeplot(pll, Float64.(f[:]))
function pllbodeplot(pll::PLL, f::Vector{Float64})

	OL = pll["Open-loop"]
	mdB,ϕ = bode(OL,f)
	fₘᵢₙ = minimum(f)
	fₘₐₓ = maximum(f)
	mdBₘᵢₙ = minimum(mdB)
	mdBₘₐₓ = maximum(mdB)
	ϕₘᵢₙ = minimum(ϕ)
	ϕₘₐₓ = maximum(ϕ)
	rₘ = mdBₘₐₓ-mdBₘᵢₙ
	rᵩ = ϕₘₐₓ-ϕₘᵢₙ

	
	fig, ax = subplots(2, 1, sharex=true)
	fig.canvas.manager.set_window_title("PLL Bode plot")
	ax[1].semilogx(f, mdB, lw=2)
	ax[2].semilogx(f, ϕ, lw=2)
	ax[1].set_ylabel("Amplitude (dB)")
	ax[2].set_ylabel("Phase (degrees)")
	ax[1].grid(b=true, which="major", axis="both")
	ax[2].grid(b=true, which="major", axis="both")
	ax[1].grid(b=true, which="minor", axis="both", linestyle=":")
	ax[2].grid(b=true, which="minor", axis="both", linestyle=":")
	ax[2].set_xlabel("Frequency (Hz)")
	ax[1].set_xlim(fₘᵢₙ, fₘₐₓ)
	ax[2].set_xlim(fₘᵢₙ, fₘₐₓ)
	ax[1].set_ylim(mdBₘᵢₙ-0.1*rₘ, mdBₘₐₓ+0.1*rₘ)
	ax[2].set_ylim(ϕₘᵢₙ-0.1*rᵩ, ϕₘₐₓ+0.1*rᵩ)

	try
		info = pllbodeinfo(pll)
		f₀ = info.ϕₘ.first
		pₘ = info.ϕₘ.second
		f₀s,f₀p = num2si(f₀)
		pₘs,pₘp = num2si(pₘ)
		ax[1].text(f₀,0,@sprintf("%.1f%sHz",f₀s,f₀p), size=8, ha=:left, va=:bottom)
		ax[2].text(f₀,pₘ-180,@sprintf("%.1f%s°",pₘs,pₘp), size=8, ha=:left, va=:bottom)
		ax[2].annotate("", xytext=(f₀,-180), xy=(f₀, pₘ-180), arrowprops=Dict("arrowstyle"=>"->"))
		ax[1].axhline([0],    c="black", lw=1)
		ax[2].axhline([-180], c="black", lw=1)
		ax[1].axvline([f₀], linestyle="--")
	catch
	end
	ax[1].xaxis.set_major_formatter(matplotlib.ticker.EngFormatter())
	ax[2].xaxis.set_major_formatter(matplotlib.ticker.EngFormatter())
	
	display(gcf())
	return
end