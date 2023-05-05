
using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")
elements, nodes = import_quad("./msh/patch_test_quad.msh")

Ē = 1.0
ν̄ = 0.49999999
# ν̄ = 0.3

E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵛ"])
set∇𝝭!(elements["Ωᵛ"])
set𝝭!(elements["Γᵍ"])

# u(x,y) = x^2+y^2
# ∂u∂x(x,y) = 2*x
# ∂u∂y(x,y) = 2*y
# v(x,y) = x^2-2*x*y
# ∂v∂x(x,y) = 2*x-2*y
# ∂v∂y(x,y) = -2*x

u(x,y) = x+y
∂u∂x(x,y) = 1.0
∂u∂y(x,y) = 1.0
v(x,y) = x+y
∂v∂x(x,y) = 1.0
∂v∂y(x,y) = 1.0

# u(x,y) = x*y
# ∂u∂x(x,y) = y
# ∂u∂y(x,y) = x
# v(x,y) = x*y
# ∂v∂x(x,y) = y
# ∂v∂y(x,y) = x

prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

op_Ω = Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν)
op_Γ = Operator{:∫vᵢgᵢds}(:α=>1e15*E)

op_Ωᵛ = Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄)
op_Ωᵈ = Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄)

k = zeros(2*nₚ,2*nₚ)
kᵛ = zeros(2*nₚ,2*nₚ)
kᵈ = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

op_Ω(elements["Ω"],k)
op_Γ(elements["Γᵍ"],kα,f)
op_Ωᵛ(elements["Ωᵛ"],kᵛ)
op_Ωᵈ(elements["Ω"],kᵈ)

# d = (k+kα)\f
d = (kᵛ+kᵈ+kα)\f
# d = (kᵛ+kα)\f

d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)

op = Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ω"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
h1,l2 = op(elements["Ω"])

rk = rank(k)
rkᵛ = rank(kᵛ)
rkᵈ = rank(kᵈ)