# ------------------------------------------------------------------------------------------
# Imported types

using AbstractTrees




#-------------------------------------------------------------------------------------------
# Internal functions

export dB10, dB20
dB10(x) = 10*log10(x)
dB20(x) = 2*dB10(x)




#-------------------------------------------------------------------------------------------
# PLL type & constructors
export PLL
"""
	PLL

Linear time invariant (LTI) PLL model, based on a continuous-time approximation [1]. The PLL
model has the following structure:

       ┏━━━━━━┓   ┏━━━┓   ┏━━━━━━━━┓       ┏━━━━━━┓
    ───┨ I(s) ┠───┨ + ┠───┨  L(s)  ┠───┬───┨ O(s) ┠───
       ┗━━━━━━┛   ┗━┯━┛   ┗━━━━━━━━┛   │   ┗━━━━━━┛
                  ─ │     ┏━━━━━━━━┓   │
                    └─────┨  F(s)  ┠───┘
                          ┗━━━━━━━━┛

where I(s), L(s), F(s) and O(s) represent the input, feedforward, feedback and output 
transfer functions, respectively. An additional closed-loop transfer function, G(s), is 
formed by the feedback network of L(s) and F(s):

       L(s)
    ───────────
    1+L(s)⋅F(s)

Each branch is made up of a set of zero or more instances of [`Block`](@ref). 

[1] TODO! Get reference
"""
struct PLL
	name::String
	description::String
	input::Vector{CascadedBlock}
	forward::Vector{CascadedBlock}
	feedback::Vector{CascadedBlock}
	output::Vector{CascadedBlock}

	_ol_tf::TF
	_cl_tf::TF
	_in_tf::TF
	_out_tf::TF
	_fwd_tf::TF
	_fbk_tf::TF

	@doc """
		PLL(name, input, forward, feedback, output, description="")

	Create a PLL model.

	The input, forward, feedback and output blocks are given by the respective iterable.
	Note that ordered iterables must be used to guarantee the block order. 

	# Arguments
	- `name::String`: model name
	- `input`: ordered iterable of blocks on the input path
	- `forward`: ordered iterable of blocks on the forward path
	- `feedback`: ordered iterable of blocks on the feedback path
	- `output`: ordered iterable of blocks on the output path
	- `description::String`: model description
	"""
	function PLL(name::String, input, forward, feedback, output, description::String="")
		hin  = reduce(*, (s->s.tf).(input),     init=TF())
		hfwd = reduce(*, (s->s.tf).(forward),   init=TF())
		hfbk = reduce(*, (s->s.tf).(feedback),  init=TF())
		hout = reduce(*, (s->s.tf).(output),    init=TF())
		hol = hfwd*hfbk
		hcl = hfwd/(1+hol)
		inp = _parse_branch(input, hcl*hout)
		fwd = _parse_branch(forward, hout/(1+hol))
		fbk = _parse_branch(feedback, -hout*hcl)
		out = _parse_branch(output)
		new(name, description, inp, fwd, fbk, out, hol, hcl, hin, hout, hfwd, hfbk)
	end
end


function _parse_branch(blocks, tf::TF=TF())
	vec = Vector{CascadedBlock}(undef, length(blocks))
	ntf = tf
	nᵢ = length(vec)
	for bᵢ in reverse(eachindex(blocks))
		vec[nᵢ] = CascadedBlock(blocks[bᵢ], TF(ntf,name=blocks[bᵢ].name*" NTF"))
		nᵢ  -= 1
		ntf *= blocks[bᵢ].tf
	end
	return vec
end


"""
	PLL()

Create an example PLL model.
"""
function PLL()
	fref = 80e6
	Nref = 8
	Icp = 1e-3
	Kvco = 10e6
	Nmmd = 2e9*Nref/fref
	Nout = 2
	fnoise = 10 .^(1:8)
	pnoise = Dict(
		"ref"  => -[120,130,140,150,160,165,166,166],
		"rdiv" => -[130,140,150,160,164,165,165,165],
		"vco"  => -[000,030,060,080,100,120,140,160]
	)
	ref  = pllref("REF", PinkNoise("Ref. clock",fnoise, pnoise["ref"]), "80MHz reference clock")
	rdiv = plldiv("RDIV", Nref, PinkNoise("Ref. div.", fnoise, pnoise["ref"]),"Reference divider")
	pfd = pllgain("PFD", 1/(2*pi))
	chp = pllgain("CP", Icp)
	# lpf = pllfilter("LPF", [9225,8547], [19,540,19]*1e-12)
	lpf, = pllfilter("LPF", Icp*Kvco/Nmmd, 100e3, 53, 300, 0.2)
	vco = pllvco("VCO", Kvco, PinkNoise("VCO", fnoise, pnoise["vco"]))
	mmd = pllmmd("MMD", Nmmd, fref/Nref, 3)
	odiv = plldiv("ODIV", Nout)

	input = [ref,rdiv]
	forward = [pfd,chp,lpf,vco]
	feedback = [mmd]
	output = [odiv]
	pll = PLL("Example PLL model", input, forward, feedback, output, "Generated by the PllToolbox module")
	# fref = 1e6
	# Nref = 1
	# Icp = 1e-3
	# Kvco = 1e6
	# Nmmd = 103
	# Nout = 1
	# ref = pllref("80MHz reference clock", PinkNoise("Reference noise",[100,1e6,1e8,1e9],[-20,-130.0,-160.0,-170.0]))
	# rdiv = plldiv("Reference divider", Nref, WhiteNoise("Reference divider noise", 10^-13))
	# pfd = pllgain("PFD", 1/(2*pi), WhiteNoise("PFD noise", 10^-12))
	# chp = pllgain("Charge-pump", Icp, WhiteNoise("Charge-pump noise", 10^-15))
	# lpf, = pllfilter("LPF", [1e3], [1.0,100.0]*1e-9)
	# vco = pllvco("2GHz VCO", Kvco, PinkNoise("VCO noise",[100,1e6,1e8,1e9],[0.0,-110.0,-140.0,-145.0]))
	# mmd = pllmmd("MMD", Nmmd, fref/Nref, 3)
	# odiv = plldiv("Output divider", Nout, WhiteNoise("Output divider noise", 10^-13))
	# pll = PLL("Example PLL model", "Generated by the PllToolbox module")
	# addinput!(pll, ref, rdiv)
	# addforward!(pll, pfd, chp, lpf, vco)
	# addfeedback!(pll, mmd)
	# addoutput!(pll, odiv)
	return pll
end




#-------------------------------------------------------------------------------------------
# PLL transfer functions
"""
    getindex(pll::PLL, key::AbstractString)

Return the PLL transfer function identified by `key` (case insensitive).
Throws a `KeyError` unless `key` is one of:

	"closed-loop", "cl"
	"open-loop", "ol"
	"input", "in"
	"output", "out"
	"forward", "fwd"
	"feedback", "fbk"

"""
function getindex end

Base.getindex(pll::PLL, keys::AbstractString...) = getindex.(Ref(pll),keys)

function Base.getindex(pll::PLL, key::AbstractString)
	lowerkey = lowercase(key)
	clkeys  = ("cl","closed-loop")
	olkeys  = ("ol","open-loop")
	inkeys  = ("in","input")
	outkeys = ("out","output")
	fwdkeys = ("fwd","forward")
	fbkkeys = ("fbk","feedback")
	if any(lowerkey.==olkeys)   return getfield(pll, :_ol_tf)   end
	if any(lowerkey.==clkeys)   return getfield(pll, :_cl_tf)   end
	if any(lowerkey.==inkeys)   return getfield(pll, :_in_tf)   end
	if any(lowerkey.==outkeys)  return getfield(pll, :_out_tf)  end
	if any(lowerkey.==fwdkeys)  return getfield(pll, :_fwd_tf)  end
	if any(lowerkey.==fbkkeys)  return getfield(pll, :_fbk_tf)  end
	throw(KeyError(key))
end




export pllntf
"""
	ntf = pllntf(pll)

Return an array of noise transfer functions (NTFs) for a PLL. 
"""
function pllntf(pll::PLL)
	[cb.ntf for cb=vcat(pll.input,pll.forward,pll.feedback,pll.output)]
end




#-------------------------------------------------------------------------------------------
# Interface implementations
function Base.show(io::IO, ::MIME"text/plain", x::PLL)
	println(io, "PLL structure")
	println(io, "    name: $(x.name)")
	if length(x.input)>0
		println(io, "    input blocks:")
		for b in x.input
			println(io, "        $(b.name)")
		end
	else
		println(io, "    input blocks: none")
	end
	if length(x.forward)>0
		println(io, "    forward blocks:")
		for b in x.forward
			println(io, "        $(b.name)")
		end
	else
		println(io, "    forward blocks: none")
	end
	if length(x.feedback)>0
		println(io, "    feedback blocks:")
		for b in x.feedback
			println(io, "        $(b.name)")
		end
	else
		println(io, "    feedback blocks: none")
	end
	if length(x.output)>0
		println(io, "    output blocks:")
		for b in x.output
			println(io, "        $(b.name)")
		end
	else
		println(io, "    output blocks: none")
	end
end



struct PathNode
	name::String
	children::Vector{CascadedBlock}
end


function AbstractTrees.children(node::PLL)
	(
		PathNode("Input", node.input),
		PathNode("Forward", node.forward),
		PathNode("Feedback", node.feedback),
		PathNode("Output", node.output)
	)
end
AbstractTrees.children(node::PathNode) = node.children


AbstractTrees.printnode(io::IO, node::PLL) = print(io, node.name)
AbstractTrees.printnode(io::IO, node::PathNode) = print(io, node.name)
