
using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")
fid_𝑢 = "./msh/patch_test_10.msh"
fid_𝑝 = "./msh/patch_test_10.msh"
elements, nodes = import_rkgsi_fem(fid_𝑢,fid_𝑝)
# elements, nodes = import_rkgsi("./msh/patch_test_10.msh")

Ē = 1.0
# ν̄ = 0.49999999
ν̄ = 0.3

E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

nₚ = length(nodes)
nₑ = length(elements["Ω"])

s = 2.5/10*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω̃"])
set𝝭!(elements["Ωᶠ"])
set∇𝝭!(elements["Ωᶠ"])
set∇𝝭!(elements["Ω̄"])
set∇𝝭!(elements["Γᵍ"])
set∇𝝭!(elements["Ωᵉ"])

# 𝗚 = elements["Ω̄"][1].𝗚
# 𝗴₁ = elements["Ω̄"][1].𝗴₁
# 𝗴₂ = elements["Ω̄"][1].𝗴₂
# m = zeros(nₚ,nₚ)
# op = Operator{:∫vudΩ}()
# op(elements["Ωᶠ"],m)
# err = m-𝗚
# err0 = 0.0
# for a in elements["Ω̄"]
#     𝓒 = a.𝓒
#     𝓖 = a.𝓖
#     for ξ in 𝓖
#         𝑤 = ξ.𝑤
#         B₁ = ξ[:∂𝝭∂x]
#         B₂ = ξ[:∂𝝭∂y]
#         temp = 0.0
#         for (i,xᵢ) in enumerate(𝓒)
#             temp += B₁[i]*xᵢ.x^0
#         end
#         println(temp)
#         global err0 += temp^2*𝑤
#     end
# end

u(x,y) = x+y
∂u∂x(x,y) = 1.0
∂u∂y(x,y) = 1.0
v(x,y) = x+y
∂v∂x(x,y) = 1.0
∂v∂y(x,y) = 1.0

prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

op_Ω = Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν)
op_Γ = Operator{:∫σᵢⱼnⱼgᵢds}(:E=>E,:ν=>ν)

op_Ωᵛ = Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄)
op_Ωᵈ = Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄)

k = zeros(2*nₚ,2*nₚ)
kᵛ = zeros(2*nₚ,2*nₚ)
kᵈ = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

op_Ω(elements["Ω̃"],k)
op_Γ(elements["Γᵍ"],kα,f)
op_Ωᵛ(elements["Ω̄"],kᵛ)
op_Ωᵈ(elements["Ω̃"],kᵈ)

# d = (k+kα)\f
d = (kᵛ+kᵈ+kα)\f
# d = (kᵛ+kα)\f

d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)

op = Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
prescribe!(elements["Ωᵉ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵉ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵉ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵉ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵉ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵉ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
h1,l2 = op(elements["Ωᵉ"])

rk = rank(k)
rkᵛ = rank(kᵛ)
rkᵈ = rank(kᵈ)
