using Revise, ApproxOperator, TOML, LinearAlgebra

config = TOML.parsefile("./toml/rk13.toml")
file1 = "./msh/square_bubble_24.msh"
file2 = "./msh/heat_12.msh"
# elements, nodes = importmsh("./msh/square_rk6_1.msh",config)
# elements, nodes = importmsh(file1,config)
elements, nodes = importmsh(file1,file2,config)

dgt = 10
dofs = Tuple{Float64,Float64}[]
for elm in elements["Ω"]
    𝓖 = elm.𝓖
    for ξ in 𝓖
        x = trunc(ξ.x; digits=dgt)
        y = trunc(ξ.y; digits=dgt)
        push!(dofs,(x,y))
    end
end
unique!(dofs)


nᵢ = length(dofs)
k = zeros(nᵢ,nᵢ)
g₁ᵢ = zeros(13,6)
g₂ᵢ = zeros(13,6)
g₁ⱼ = zeros(13,6)
g₂ⱼ = zeros(13,6)
G⁻¹ = [  36. -120. -120.   90.  180.   90.;
       -120.  600.  300. -540. -720. -180.;
       -120.  300.  600. -180. -720. -540.;
         90. -540. -180.  540.  540.   90.;
        180. -720. -720.  540. 1440.  540.;
         90. -180. -540.   90.  540.  540.]
get𝒒(ξ,η) = (1.0,ξ,η,ξ^2,ξ*η,η^2)
get∂𝒒∂ξ(ξ,η) = (0.0,1.0,0.0,2*ξ,η,0.0)
get∂𝒒∂η(ξ,η) = (0.0,0.0,1.0,0.0,ξ,2*η)

for elm in elements["Ω"][1:end]
    𝓖 = elm.𝓖
    for (i,ξᵢ) in enumerate(𝓖)
        ξ = ξᵢ.ξ
        η = ξᵢ.η
        q = get𝒒(ξ,η)
        ∂q∂ξ = get∂𝒒∂ξ(ξ,η)
        ∂q∂η = get∂𝒒∂η(ξ,η)
        wb = ξᵢ.wᵇ
        wi = ξᵢ.w
        D₁ = ξᵢ.D₁
        D₂ = ξᵢ.D₂
        D₁₁ = ξᵢ.D₁₁
        D₂₁ = ξᵢ.D₂₁
        D₁₂ = ξᵢ.D₁₂
        D₂₂ = ξᵢ.D₂₂
        for (k,qₖ) in enumerate(q)
            ∂qₖ∂ξ = ∂q∂ξ[k]
            ∂qₖ∂η = ∂q∂η[k]
            g₁ᵢ[i,k] = wi/2.0*(∂qₖ∂ξ*D₁₁+∂qₖ∂η*D₂₁)+wb*qₖ*D₁
            g₂ᵢ[i,k] = wi/2.0*(∂qₖ∂ξ*D₁₂+∂qₖ∂η*D₂₂)+wb*qₖ*D₂
        end
    end
    for (j,ξⱼ) in enumerate(𝓖)
        ξ = ξⱼ.ξ
        η = ξⱼ.η
        q = get𝒒(ξ,η)
        ∂q∂ξ = get∂𝒒∂ξ(ξ,η)
        ∂q∂η = get∂𝒒∂η(ξ,η)
        wb = ξⱼ.wᵇ
        wi = ξⱼ.w
        D₁ = ξⱼ.D₁
        D₂ = ξⱼ.D₂
        D₁₁ = ξⱼ.D₁₁
        D₂₁ = ξⱼ.D₂₁
        D₁₂ = ξⱼ.D₁₂
        D₂₂ = ξⱼ.D₂₂
        for (k,qₖ) in enumerate(q)
            ∂qₖ∂ξ = ∂q∂ξ[k]
            ∂qₖ∂η = ∂q∂η[k]
            g₁ⱼ[j,k] = wi/2.0*(∂qₖ∂ξ*D₁₁+∂qₖ∂η*D₂₁)+wb*qₖ*D₁
            g₂ⱼ[j,k] = wi/2.0*(∂qₖ∂ξ*D₁₂+∂qₖ∂η*D₂₂)+wb*qₖ*D₂
        end
    end
    𝐴 = 𝓖[1].𝐴
    for (i,ξᵢ) in enumerate(𝓖)
        x = trunc(ξᵢ.x; digits=dgt)
        y = trunc(ξᵢ.y; digits=dgt)
        I = findfirst(x_->x_==(x,y),dofs)
        for (j,ξⱼ) in enumerate(𝓖)
            x = trunc(ξⱼ.x; digits=dgt)
            y = trunc(ξⱼ.y; digits=dgt)
            J = findfirst(x_->x_==(x,y),dofs)
            for ii in 1:6
                for jj in 1:6
                    k[I,J] += g₁ᵢ[i,ii]*G⁻¹[ii,jj]*g₁ⱼ[j,jj]/𝐴 + g₂ᵢ[i,ii]*G⁻¹[ii,jj]*g₂ⱼ[j,jj]/𝐴
                    # k[I,J] += g₁ᵢ[i,ii]*G⁻¹[ii,jj]*g₁ⱼ[j,jj]/𝐴
                end
            end
        end
    end
end

rk = rank(k)

nₚ = length(nodes)
s = 3.5/12*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
set𝝭!(elements["Ω"])
Ψ = zeros(nₚ,nᵢ)
for elm in elements["Ω"]
    𝓒 = elm.𝓒
    𝓖 = elm.𝓖
    for (j,ξ) in enumerate(𝓖)
        N = ξ[:𝝭]
        x = trunc(ξ.x; digits=dgt)
        y = trunc(ξ.y; digits=dgt)
        J = findfirst(x_->x_==(x,y),dofs)
        for (i,x) in enumerate(𝓒)
            I = x.𝐼
            Ψ[I,J] = N[i]
        end
    end
end

rs = rank(Ψ)

rsks = rank(Ψ*k*Ψ')
println("rank of G⁻¹ = $rk, rank of Ψ = $rs, rank of ΨG⁻¹ψᵀ = $rsks")