
using Revise,ApproxOperator,LinearAlgebra
include("input.jl")

fid_𝑢 = "./msh/cantilever_8.msh"
fid_𝑝 = "./msh/cantilever_8.msh"
elements, nodes, nodes_𝑝 = import_rkgsi_mix_quadratic(fid_𝑢,fid_𝑝)

nᵤ = length(nodes)
nₚ = length(nodes)
sᵤ = 2.5*12/8*ones(nᵤ)
push!(nodes,:s₁=>sᵤ,:s₂=>sᵤ,:s₃=>sᵤ)
sₚ = 2.5*12/8*ones(nₚ)
push!(nodes_𝑝,:s₁=>sₚ,:s₂=>sₚ,:s₃=>sₚ)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω̃"])
set𝝭!(elements["Ωˢᵖ"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Ω̃ᵖ"])
set∇𝝭!(elements["Ω̄"])
set𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Γᵍ"])

# temp = 0.0
# 𝗚 = elements["Ω̄"][1].𝗚
# 𝗴₁ = elements["Ω̄"][1].𝗴₁
# 𝗴₂ = elements["Ω̄"][1].𝗴₂
# for (i,pᵢ) in enumerate(nodes_𝑝)
#     # for (j,pⱼ) in enumerate(nodes_𝑝)
#     for (j,pⱼ) in enumerate(nodes)
#         # global temp += pᵢ.x*𝗚[i,j]*pⱼ.x
#         global temp += pᵢ.x*𝗴₂[i,j]*pⱼ.y
#     end
# end
# temp-13824
# LinearAlgebra.cond(𝗚)
# for a in elements["Ωˢᵖ"][[1]]
# for a in elements["Ωᵖ"]
# for a in elements["Ω̃ᵖ"]
for a in elements["Ω̄"]
    for ξ in a.𝓖
        # 𝝭 = ξ[:𝝭]
        ∂𝝭∂x = ξ[:∂𝝭∂x]
        ∂𝝭∂y = ξ[:∂𝝭∂y]
        u = 0.0
        for (i,x) in enumerate(a.𝓒)
            # u += 𝝭[i]*x.y
            u += ∂𝝭∂x[i]
        end
        # if abs(u - ξ.y) > 1e-13
        if abs(u - 0.0) > 1e-13
        # if abs(u) > 1e-13
            println(u)
            println(1.0)
            error("consistency condition is not satisfied")
        end
    end
end
# for (a,b) in zip(elements["Ωˢᵖ"],elements["Ωᵖ"])
#     𝓒 = a.𝓒
#     for (ξa,ξb) in zip(a.𝓖,b.𝓖)
#         𝝭a = ξa[:𝝭]
#         𝝭b = ξb[:𝝭]
#         for (i,xᵢ) in enumerate(𝓒)
#             if 𝝭a[i] ≠ 𝝭b[i]
#                 println(𝝭a[i])
#                 println(𝝭b[i])
#                 error("shape function is not equal")
#             end
#         end
#     end
# end

# E = 3E6;ν = 0.3;P = 1000;L = 48;D = 12;
# I = D^3/12
# EI = E*I

# prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.)
# prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
# prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
# prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
# prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
# prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
# prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

# coefficient = (:E=>E,:ν=>ν,:α=>1E3*E)
# ops = [
#     Operator{:∫∫εᵢⱼσᵢⱼdxdy}(coefficient...),
#     Operator{:∫vᵢtᵢds}(coefficient...),
#     Operator{:∫σᵢⱼnⱼgᵢds}(coefficient...),
#     Operator{:∫vᵢgᵢds}(coefficient...),
#     Operator{:Hₑ_PlaneStress}(coefficient...)
# ]

# k = zeros(2*nₚ,2*nₚ)
# f = zeros(2*nₚ)

# ops[1](elements["Ω̃"],k)
# # ops[1](elements["Ω̄"],k)
# ops[2](elements["Γᵗ"],f)
# ops[3](elements["Γᵍ"],k,f)
# ops[4](elements["Γᵍ"],k,f) 

# d = k\f
# d₁ = d[1:2:2*nₚ]
# d₂ = d[2:2:2*nₚ]
# push!(nodes,:d₁=>d₁,:d₂=>d₂)

# set∇𝝭!(elements["Ωₑ"])
# prescribe!(elements["Ωₑ"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
# prescribe!(elements["Ωₑ"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
# prescribe!(elements["Ωₑ"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
# prescribe!(elements["Ωₑ"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
# prescribe!(elements["Ωₑ"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
# prescribe!(elements["Ωₑ"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)
# h1,l2 = ops[5](elements["Ωₑ"])
