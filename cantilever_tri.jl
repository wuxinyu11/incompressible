
using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")
elements, nodes = import_tri3("./msh/cantilever_8.msh")

Ē = 3e6
# ν̄ = 0.499999999999999
ν̄ = 0.3

E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
# E = 3e6
# ν = 0.3

nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γᵍ"])

P = 1000;L = 48;D = 12;
I = D^3/12
EI = E*I
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

op_Ω = Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν)
op_Γᵍ = Operator{:∫vᵢgᵢds}(:α=>1e10*E)
op_Γᵗ = Operator{:∫vᵢtᵢds}()

k = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
fα = zeros(2*nₚ)

op_Ω(elements["Ω"],k)
op_Γᵗ(elements["Γᵗ"],f)
op_Γᵍ(elements["Γᵍ"],kα,fα)

d = (k+kα)\(f+fα)

d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)

set𝝭!(elements["Ωᵉ"])
set∇𝝭!(elements["Ωᵉ"])
op = Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
prescribe!(elements["Ωᵉ"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Ωᵉ"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Ωᵉ"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
prescribe!(elements["Ωᵉ"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
prescribe!(elements["Ωᵉ"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
prescribe!(elements["Ωᵉ"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)
h1,l2 = op(elements["Ωᵉ"])