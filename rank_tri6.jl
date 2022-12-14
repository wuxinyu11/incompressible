
using Revise, ApproxOperator, TOML, LinearAlgebra

config = TOML.parsefile("./toml/tri6.toml")
elements, nodes = importmsh("./msh/square_tri6_1.msh",config)

nₚ = length(nodes)

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ"])

prescribe!(elements["Γ"],:g=>(x,y,z)->1.0+2.0*x+3.0*y)
op_Ω = Operator{:∫∇v∇udΩ}(:k=>1.0)
op_Γ = Operator{:∫vgdΓ}(:α=>1e8)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
op_Ω(elements["Ω"][1:1],k)
# op_Γ(elements["Γ"],k,f)

rk = rank(k)
