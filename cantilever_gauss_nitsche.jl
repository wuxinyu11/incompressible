
using  ApproxOperator, LinearAlgebra, Printf

include("input.jl")

for n in 2:5 
println(n-1)
    ndiv = 2^n
#  ndiv = 32

elements, nodes = import_gauss_quadratic("./msh/cantilever_"*string(ndiv)*".msh",:TriGI13)
nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 2.5*12.0/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Γᵍ"])
set∇𝝭!(elements["Ωᵉ"])
P = 1000
Ē = 3e6
ν̄ = 0.499999999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
L = 48
D = 12
I = D^3/12
EI = E*I

prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3*x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
# op_Ω = Operator(:∫∫εᵢⱼσᵢⱼdxdy,:E=>E,:ν=>ν)
# op_Ωᵛ = Operator(:∫∫εᵛᵢⱼσᵛᵢⱼdxdy,:E=>Ē,:ν=>ν̄)
# op_Ωᵈ = Operator(:∫∫εᵈᵢⱼσᵈᵢⱼdxdy,:E=>Ē,:ν=>ν̄)
# op_Γᵗ = Operator(:∫vᵢtᵢds)
# op_Γᵍ = Operator(:∫vᵢgᵢds,:α=>1e7*Ē)
# op_He = Operator(:Hₑ_PlaneStress,:E=>E,:ν=>ν,:α=>1e7*Ē)

# ops = [Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν)
#        Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄)
#        Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄)
#        Operator{:∫vᵢtᵢds}()
#        Operator{:∫vᵢgᵢds}(:α=>1e7*Ē)
#        Operator{:Hₑ_PlaneStress,}(:E=>E,:ν=>ν)]
# ops = [
    
#     Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2}(:E=>Ē,:ν=>ν̄),
#     Operator{:∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2}(:E=>Ē,:ν=>ν̄),
#     Operator{:∫vᵢtᵢds}(),
#     Operator{:∫σᵢⱼnⱼgᵢds}(:E=>E,:ν=>ν),
#     Operator{:∫vᵢgᵢds}(:α=>1e15*E),
#     Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
# ]

# coefficient = (:E=>E,:ν=>ν,:α=>1E7*E)
ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
       Operator{:∫vᵢtᵢds}(),
    #    Operator{:∫σᵢⱼnⱼgᵢds}(:E=>E,:ν=>ν),
       Operator{:∫vᵢgᵢds}(:α=>1e3*E),
       Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
]

kα = zeros(2*nₚ,2*nₚ)
k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
fα = zeros(2*nₚ)
d = zeros(2*nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)

push!(nodes,:d₁=>d₁,:d₂=>d₂)

ops[1](elements["Ω"],k)
ops[2](elements["Γᵗ"],f)
ops[3](elements["Γᵍ"],k,f)
# ops[3](elements["Γᵍ"],kα,fα)
# ops[4](elements["Γᵍ"],k,f)
d = k\f
# d = (k+kα)\(f+fα)
d₁ .= d[1:2:2*nₚ]
d₂ .= d[2:2:2*nₚ]

prescribe!(elements["Ωᵉ"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Ωᵉ"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Ωᵉ"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
prescribe!(elements["Ωᵉ"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
prescribe!(elements["Ωᵉ"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
prescribe!(elements["Ωᵉ"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)
# h1,l2 = ops[5](elements["Ωᵉ"])
h1,l2 = ops[4](elements["Ωᵉ"])
println(h1)
println(l2)
end
        