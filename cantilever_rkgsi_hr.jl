
using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")

ndiv = 8

fid_𝑢 = "./msh/cantilever_"*string(ndiv)*".msh"
fid_𝑝 = "./msh/cantilever_"*string(ndiv)*".msh"

elements, nodes = import_rkgsi_fem(fid_𝑢,fid_𝑝)

nₚ = length(nodes)
nₑ = length(elements["Ω"])

#  s = 3*12 / ndiv * ones(nₚ)
 s = 2.5*12 / ndiv * ones(nₚ)


 push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω̃"])
set𝝭!(elements["Ωᶠ"])
set∇𝝭!(elements["Ωᶠ"])
set∇𝝭!(elements["Ω̄"])
set∇𝝭!(elements["Γᵍ"])
set∇𝝭!(elements["Ωᵉ"])

P = 1000
Ē = 3e6
# ν̄ = 0.499999999
ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
L = 48
D = 12
I = D^3/12
EI = E*I
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

coefficient = (:E=>E,:ν=>ν)
ops = [Operator{:∫∫εᵢⱼσᵢⱼdxdy}(coefficient...),
       Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄),
       Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄),
       Operator{:∫vᵢtᵢds}(coefficient...),
       Operator{:∫σᵢⱼnⱼgᵢds}(coefficient...),
       Operator{:∫σ̄ᵢⱼnⱼgᵢds}(coefficient...),
       Operator{:Hₑ_PlaneStress}(coefficient...)]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

# ops[1](elements["Ω̃"],k)
ops[3](elements["Ω̃"],k)
ops[2](elements["Ω₁"],k)
ops[4](elements["Γᵗ"],f)
ops[5](elements["Γᵍ"],k,f)
ops[6](elements["Γᵍ"],k,f)

 d = k\f

d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
set𝓖!(elements["Ω"],:TriGI16,:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)

# set𝓖!(elements["Ω"],:TriGI16,:∂1,:∂x,:∂y,:∂z)
# set𝝭!(elements["Ω"])
set∇₂𝝭!(elements["Ω"])
prescribe!(elements["Ω"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Ω"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)
h1,l2 = ops[7](elements["Ω"])
# L2 = log10(l2)
# h1 = log10(h1)
# h = log10(12.0/ndiv)
 L2 = log10(l2)

 index = [8,16,32,64]
XLSX.openxlsx("./xlsx/cantilever.xlsx", mode="rw") do xf
    row = "G"
    𝐿₂ = xf[2]
    #  𝐻₁ = xf[3]
#     𝐻₂ = xf[4]
#     𝐻₃ = xf[5]
    ind = findfirst(n->n==ndiv,index)+1
    row = row*string(ind)
    𝐿₂[row] = log10(l2)
    #  𝐻₁[row] = log10(h1)
#     𝐻₂[row] = log10(h2)
#     𝐻₃[row] = log10(h3)
end