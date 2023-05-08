
using  ApproxOperator, LinearAlgebra, Printf
include("input.jl")
elements, nodes = import_gauss_quadratic("./msh/cook_membrance_10.msh",:TriGI3)

κ = 400942
μ = 80.1938
E = 9*κ*μ/(3*κ+μ)
ν = (3*κ-2*μ)/2/(3*κ+μ)
# E = 70.0
# ν = 0.3333

nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = 2.5*44/10*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γᵍ"])

prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->100.0)

prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫σᵢⱼnⱼgᵢds}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e4*E)
]

k = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
fα = zeros(2*nₚ)
fint = zeros(2*nₚ)
fext = zeros(2*nₚ)
d = zeros(2*nₚ)
Δd= zeros(2*nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)

push!(nodes,:d₁=>d₁,:d₂=>d₂)
ops[1](elements["Ω"],k)
ops[2](elements["Γᵗ"],f)
ops[3](elements["Γᵍ"],k,f)
ops[4](elements["Γᵍ"],k,f)
d = k\f
        d₁ .= d[1:2:2*nₚ]
        d₂ .= d[2:2:2*nₚ]

u₁=d₁[3]
u₂=d₂[3]
println(u₂)
