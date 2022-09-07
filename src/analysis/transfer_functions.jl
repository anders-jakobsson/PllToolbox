export plltf
"""
	G,L,F,I,O = plltf(pll)

Return the total, forward, feedback, input and output transfer functions of a PLL.
"""
function plltf(pll::PLL)
	I = reduce(*, (s->s.H).(pll.input),    init=tf(1))
	L = reduce(*, (s->s.H).(pll.forward),  init=tf(1))
	F = reduce(*, (s->s.H).(pll.feedback), init=tf(1))
	O = reduce(*, (s->s.H).(pll.output),   init=tf(1))
	G = L/(1+L*F)
	(G,L,F,I,O)
end




export pllntf
"""
	ntf = pllntf(pll)

Return an array of noise transfer functions (NTFs) for a PLL. 
"""
function _pllnoise_branch(blocks::Vector{Block}, f::Vector{Float64}, h::LTISystem=tf(1))
	num_noisy = count(b->!isempty(b.noise), blocks)
	bnoise = Vector{BlockNoise}(undef, num_noisy)
	htot = tf(h)
	nᵢ = length(bnoise)
	for bᵢ in reverse(eachindex(blocks))
		if !isempty(blocks[bᵢ].noise)
			bnoise[nᵢ] = BlockNoise(blocks[bᵢ], f, htot)
			nᵢ -= 1
		end
		htot = htot*blocks[bᵢ].H
	end
	bnoise
end