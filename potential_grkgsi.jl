using Revise, ApproxOperator, BenchmarkTools
include("input.jl")

fid_𝑢 = "./msh/patch_test.msh"
fid_𝑝 = "./msh/patch_test.msh"
elements, nodes, nodes_𝑝 = import_rkgsi_mix(fid_𝑢,fid_𝑝)

nᵤ = length(nodes)
nₚ = length(nodes)
sᵤ = 2.5/10*ones(nᵤ)
push!(nodes,:s₁=>sᵤ,:s₂=>sᵤ,:s₃=>sᵤ)
sₚ = 2.5/10*ones(nₚ)
push!(nodes_𝑝,:s₁=>sₚ,:s₂=>sₚ,:s₃=>sₚ)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω̃"])
set∇𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Ω̃ᵖ"])
set∇𝝭!(elements["Ω̄"])
# set𝝭!(elements["Γᵗ"])
# set∇𝝭!(elements["Γᵍ"])

# r = 1
# u(x,y) = (x+y)^r
# ∂u∂x(x,y) = r*(x+y)^abs(r-1)
# ∂u∂y(x,y) = r*(x+y)^abs(r-1)
# ∂²u∂x²(x,y) = r*(r-1)*(x+y)^abs(r-2)
# ∂²u∂x∂y(x,y) = r*(r-1)*(x+y)^abs(r-2)
# ∂²u∂y²(x,y) = r*(r-1)*(x+y)^abs(r-2)
# prescribe!(elements["Γᵍ"],:g=>(x,y,z)->u(x,y))
# prescribe!(elements["Ω"],:b=>(x,y,z)->-∂²u∂x²(x,y)-∂²u∂y²(x,y))

# ops = [
#     Operator{:∫∫∇v∇udxdy}(:k=>1.0),
#     Operator{:∫vbdΩ}(),
#     Operator{:∫vtdΓ}(),
#     Operator{:∫∇𝑛vgds}(:k=>1.0,:α=>1e3),
#     Operator{:H₁}()
# ]

# k = zeros(nₚ,nₚ)
# f = zeros(nₚ)
# # ops[1](elements["Ω̃"],k)
# ops[1](elements["Ω̄"],k)
# ops[2](elements["Ω"],f)
# # ops[3].(elements["Γᵗ"],f=f)
# ops[4](elements["Γᵍ"],k,f)

# d = k\f

# ApproxOperator.cal𝗠!(elements["Ω̄"])
# 𝗚 = elements["Ω̄"][1].𝗚
# 𝗴₁ = elements["Ω̄"][1].𝗴₁
# 𝗴₂ = elements["Ω̄"][1].𝗴₂
# t0 = sum([𝗚...])
# t10 = sum([𝗴₁...])
# t20 = sum([𝗴₂...])
# t1 = 0.0
# t11 = 0.0
# t21 = 0.0
# for pᵢ in nodes_𝑝
#     I = pᵢ.𝐼
#     xᵢ = pᵢ.x
#     yᵢ = pᵢ.y
#     for pⱼ in nodes_𝑝
#         J = pⱼ.𝐼
#         xⱼ = pⱼ.x
#         yⱼ = pⱼ.y
#         global t1 += yᵢ*𝗚[I,J]*yⱼ
#     end
#     for nⱼ in nodes
#         J = nⱼ.𝐼
#         xⱼ = nⱼ.x
#         yⱼ = nⱼ.y
#         global t11 += xᵢ*𝗴₁[I,J]*xⱼ^2
#         global t21 += yᵢ*𝗴₂[I,J]*yⱼ^2
#     end
# end

e1 = 0.0
for ap in elements["Ω̄"]
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        u = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            # u += B₁[i]*xᵢ.x
            u += B₁[i]
            # println(B₁[i])
            # u += 1.
            # u += B₁[i]*xᵢ.x
        end
        println(u)
        # global e1 += (u - ξ.x^0)*𝑤
        global e1 += u*𝑤
    end
end