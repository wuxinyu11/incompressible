
using ApproxOperator, LinearAlgebra, Printf
include("input.jl")

fid_𝑢 = "./msh/cook_membrance_30.msh"
fid_𝑝 = "./msh/cook_membrance_12.msh"
elements, nodes, nodes_𝑝 = import_rkgsi_mix_quadratic(fid_𝑢,fid_𝑝)

κ = 400942
μ = 80.1938
E = 9*κ*μ/(3*κ+μ)
ν = (3*κ-2*μ)/2/(3*κ+μ)
# E = 70.0
#  ν = 0.3333

nₚ = length(nodes)
n𝑝 = length(nodes_𝑝)
nₑ = length(elements["Ω"])
s = 2.5*44/30*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
s = 2.5*44/12*ones(nₚ)
push!(nodes_𝑝,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω̃"])
set𝝭!(elements["Ωˢᵖ"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Ω̃ᵖ"])
set∇𝝭!(elements["Ω̄"])
set𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Γᵍ"])

prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

ops = [
    Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    Operator{:∫∫EᵢⱼSᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢuᵢds}(:α=>1e15*E),
]
opsᵛ = [
    Operator{:Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    Operator{:∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
]
opsᵈ = [
    Operator{:Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    Operator{:∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
]

k = zeros(2*nₚ,2*nₚ)
kᵛ = zeros(2*nₚ,2*nₚ)
kᵈ = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
fα = zeros(2*nₚ)
fint = zeros(2*nₚ)
fintᵛ = zeros(2*nₚ)
fintᵈ = zeros(2*nₚ)
fext = zeros(2*nₚ)
d = zeros(2*nₚ)
Δd= zeros(2*nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)

push!(nodes,:d₁=>d₁,:d₂=>d₂)

nmax = 5
P = 0:6.25/nmax:6.25
tolerance=1.0e-10;maxiters=1000;
for (n,p) in enumerate(P)
    if n == 1
        continue
    end
    err_Δd = 1.0
    dnorm = 0.0
    # err_Δf = 1.0
    # fnorm = 0.0
    @printf "Load step=%i,p=%e \n" n p
    fill!(fext,0.0)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->p)
    ops[3](elements["Γᵗ"],fext)
    # fill!(k,0.0)
    # fill!(kα,0.0)
    # fill!(fα,0.0)
    # ops[1](elements["Ω"],k)
    # ops[4](elements["Γᵍ"],kα,fα)
    # k⁻¹ .= inv(k+kα)
    iter = 0
    while err_Δd>tolerance && iter<maxiters
        iter += 1
        fill!(kα,0.0)
        fill!(fα,0.0)
        ops[4](elements["Γᵍ"],kα,fα)

        fill!(k,0.0)
        fill!(fint,0.0)
        ops[1](elements["Ω̃"],k)
        ops[2](elements["Ω̃"],fint)

        fill!(kᵛ,0.0)
        fill!(fintᵛ,0.0)
        opsᵛ[1](elements["Ω̄"],kᵛ)
        opsᵛ[2](elements["Ω̄"],fintᵛ)
        # opsᵛ[1](elements["Ω"],kᵛ)
        # opsᵛ[2](elements["Ω"],fintᵛ)

        fill!(kᵈ,0.0)
        fill!(fintᵈ,0.0)
        opsᵈ[1](elements["Ω̃"],kᵈ)
        opsᵈ[2](elements["Ω̃"],fintᵈ)

        # if iter == 1
        #     Δd .= k⁻¹*(f+fα)
        # else
        #     Δd .= k⁻¹*f
        # end

        # f .= fext-fint
        # Δd .= (k+kα)\(f+fα)

        f .= fext-fintᵛ-fintᵈ
        Δd .= (kᵛ+kᵈ+kα)\(f+fα)

        # fnorm = norm(f)
        # fᵗnorm = fnorm+1.0
        # fᵗ .= f
        # λ = 2.0
        # while fᵗnorm ≥ fnorm && λ > tolerance
        #     # println(λ)
        #     fill!(fint,0.0)
        #     λ *= 0.5
        #     d₁ .= d[1:2:2*nₚ]+λ*Δd[1:2:2*nₚ]
        #     d₂ .= d[2:2:2*nₚ]+λ*Δd[2:2:2*nₚ]
        #     ops[2](elements["Ω"],fint)
        #     fᵗ = fext-fint
        #     fᵗnorm = norm(fᵗ)
        #     # println(fnorm)
        #     # println(fᵗnorm)
        # end
        # d .+= λ*Δd 

        d .+= Δd 
        d₁ .= d[1:2:2*nₚ]
        d₂ .= d[2:2:2*nₚ]

        Δdnorm = LinearAlgebra.norm(Δd)
        # Δdnorm = LinearAlgebra.norm(λ*Δd)
        dnorm += Δdnorm
        err_Δd = Δdnorm/dnorm
        # err_Δd = Δdnorm
        # Δfnorm = LinearAlgebra.norm(f+fα)
        # fnorm += Δfnorm
        # err_Δf = Δfnorm/fnorm

        # @printf "iter = %i, err_Δf = %e, err_Δd = %e \n" iter err_Δf err_Δd
        @printf "iter = %i, err_Δd = %e \n" iter err_Δd
    end
end 
u₁=d₁[3]
u₂=d₂[3]
println(u₂)
