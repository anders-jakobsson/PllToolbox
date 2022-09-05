using PllToolbox

pll = PLL()
f = 10 .^LinRange(1,8,100)

pllnoise(pll,f)