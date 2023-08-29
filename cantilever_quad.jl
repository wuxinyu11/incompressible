
using  ApproxOperator, LinearAlgebra, Printf, XLSX
 ndiv = 8
include("input.jl")
elements, nodes = import_quad("./msh/cantilever_quad_square_"*string(ndiv)*".msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])


set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵛ"])
set∇𝝭!(elements["Ωᵛ"])
set𝝭!(elements["Ωᵉ"])
set∇𝝭!(elements["Ωᵉ"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γᵍ"])
P = 1000
 Ē = 3e6
ν̄ = 0.49999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
L = 12
D = 12
I = D^3/12
EI = E*I
I = D^3/12
EI = E*I
# set𝒏!(elements["Γᵍ"])
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3*x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
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
       Operator{:∫vᵢgᵢds}(:α=>1e9*E),
       Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
]
opsᵛ = [
    Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]
opsᵈ = [
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]
# k = zeros(2*nₚ,2*nₚ)
kᵛ = zeros(2*nₚ,2*nₚ)
kᵛ_ = zeros(2*nₚ,2*nₚ)
kᵈ = zeros(2*nₚ,2*nₚ)
kᵍ = zeros(2*nₚ,2*nₚ)
# kα = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
# fint = zeros(2*nₚ)
d = zeros(2*nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)

push!(nodes,:d₁=>d₁,:d₂=>d₂)
        #  ops[1](elements["Ω"],k)
        # #  ops[3](elements["Ω"],k)
        # #  ops[2](elements["Ω₁"],k)
        #  ops[4](elements["Γᵗ"],f)
        #  ops[5](elements["Γᵍ"],k,f)
        # ops[1](elements["Ω"],k)
        # ops[2](elements["Γᵗ"],f)
        #   ops[3](elements["Γᵍ"],k,f)
        # ops[4](elements["Γᵍ"],k,f)       
       

        # ops[1](elements["Ω"],k)
        # ops[2](elements["Ω"],fint)
        # ops[3](elements["Γᵗ"],f)
        # ops[4](elements["Γᵍ"],k,f)
        # ops[5](elements["Γᵍ"],k,f)
        # ops[1](elements["Ω"],k)
        opsᵛ[1](elements["Ωᵛ"],kᵛ)
        opsᵛ[1](elements["Ω"],kᵛ_)
        opsᵈ[1](elements["Ωᵛ"],kᵈ)
        ops[2](elements["Γᵗ"],f)
        ops[3](elements["Γᵍ"],kᵍ,f)
        # ops[4](elements["Γᵍ"],k,f)
        # d .= (k+kα)\f
        # d = (kᵛ+kᵈ+kᵍ)\f
        d₁ .= d[1:2:2*nₚ]
        d₂ .= d[2:2:2*nₚ]
        push!(nodes,:d₁=>d₁,:d₂=>d₂)
        f = eigen(kᵈ+kᵍ,kᵛ)
        # v = eigvals(kᵈ+kᵍ,kᵛ)
        v = eigvals(kᵛ,kᵈ)
        # v = eigvals(kᵈ,kᵛ)
        # v_ = eigvals(kᵛ_,kᵈ)
        # set𝝭!(elements["Ω̄"])
        # # set∇𝝭!(elements["Ω̄"])
        #  set𝝭!(elements["Ω"])
        #  set∇𝝭!(elements["Ω"])
        # prescribe!(elements["Ω̄"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
        # prescribe!(elements["Ω̄"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
        # prescribe!(elements["Ω̄"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
        # prescribe!(elements["Ω̄"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
        # prescribe!(elements["Ω̄"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
        # prescribe!(elements["Ω̄"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)

        # prescribe!(elements["Ω"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
        # prescribe!(elements["Ω"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
        # prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
        # prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
        # prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
        # prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)
        # h1,l2 = ops[4](elements["Ω̄"])
    #     h1,l2 = ops[4](elements["Ω"])
    #    L2 = log10(l2)
    #    H1 = log10(h1)
    #    h = log10(12.0/ndiv)


# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/cantilever_quad.xlsx", mode="rw") do xf
#     Sheet = xf[4]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["B"*string(ind)] = h
#     Sheet["C"*string(ind)] = L2
#     Sheet["D"*string(ind)] = H1
# end

# d = k\f

# d₁ = d[1:2:2*nₚ]
# d₂ = d[2:2:2*nₚ]
# push!(nodes,:d₁=>d₁,:d₂=>d₂)


# set∇𝝭!(elements["Ω"])
# prescribe!(elements["Ω"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
# prescribe!(elements["Ω"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
# prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
# prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
# prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
# prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)
# h1,l2 = ops[5](elements["Ω"])
# L2 = log10(l2)
# h1 = log10(h1)
# h = log10(12.0/ndiv)
#  L2 = log10(l2)


# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/cantilever_quad.xlsx", mode="rw") do xf
    # row = "G"
    # 𝐿₂ = xf[2]
    #  𝐻₁ = xf[3]
#     𝐻₂ = xf[4]
#     𝐻₃ = xf[5]
    # ind = findfirst(n->n==ndiv,index)+1
    # row = row*string(ind)
    # 𝐿₂[row] = log10(l2)
    #  𝐻₁[row] = log10(he)
    # 𝐻₂[row] = log10(h2)
    # 𝐻₃[row] = log10(h3)
# end
