using Revise, ApproxOperator, TOML, LinearAlgebra

config = TOML.parsefile("./toml/rk6.toml")
elements, nodes = importmsh("./msh/square_rk6_1.msh",config)

nᵢ = ApproxOperator.getnᵢ(elements["Ω"])
k = zeros(nᵢ,nᵢ)
gᵢ = zeros(6,3)
gⱼ = zeros(6,3)
G⁻¹ = [9. -12. -12.;-12. 24. 12.;-12. 12. 24.]
get𝒒(ξ,η) = (1.0,ξ,η)
get∂𝒒∂ξ(ξ,η) = (0.0,1.0,0.0)
get∂𝒒∂η(ξ,η) = (0.0,0.0,1.0)
for elm in elements["Ω"][1:2]
    𝓖 = elm.𝓖
    𝐴 = ApproxOperator.get𝐴(elm)
    x₁ = elm.𝓒[1].x
    x₂ = elm.𝓒[2].x
    x₃ = elm.𝓒[3].x
    y₁ = elm.𝓒[1].y
    y₂ = elm.𝓒[2].y
    y₃ = elm.𝓒[3].y
    ∂ξ∂x = y₃-y₂
    ∂ξ∂y = x₂-x₃
    ∂η∂x = y₁-y₃
    ∂η∂y = x₃-x₁
    for (i,ξᵢ) in enumerate(𝓖)
        I = ξᵢ.𝐺
        ξ = ξᵢ.ξ
        η = ξᵢ.η
        q = get𝒒(ξ,η)
        ∂q∂ξ = get∂𝒒∂ξ(ξ,η)
        ∂q∂η = get∂𝒒∂η(ξ,η)
        wb = ξᵢ.wᵇ
        wi = ξᵢ.w
        D₁ = ξᵢ.D₁
        D₂ = ξᵢ.D₂
        for (k,qₖ) in enumerate(q)
            ∂qₖ∂ξ = ∂q∂ξ[k]
            ∂qₖ∂η = ∂q∂η[k]
            gᵢ[i,k] = wi/2.0*(∂qₖ∂ξ*∂ξ∂x+∂qₖ∂η*∂η∂x)+wb*qₖ*(D₁+D₂)
        end
    end
    for (j,ξⱼ) in enumerate(𝓖)
        J = ξⱼ.𝐺
        ξ = ξⱼ.ξ
        η = ξⱼ.η
        q = get𝒒(ξ,η)
        ∂q∂ξ = get∂𝒒∂ξ(ξ,η)
        ∂q∂η = get∂𝒒∂η(ξ,η)
        wb = ξⱼ.wᵇ
        wi = ξⱼ.w
        D₁ = ξⱼ.D₁
        D₂ = ξⱼ.D₂
        for (k,qₖ) in enumerate(q)
            ∂qₖ∂ξ = ∂q∂ξ[k]
            ∂qₖ∂η = ∂q∂η[k]
            gⱼ[j,k] = wi/2.0*(∂qₖ∂ξ*∂ξ∂x+∂qₖ∂η*∂η∂x)+wb*qₖ*(D₁+D₂)
        end
    end
    for (i,ξᵢ) in enumerate(𝓖)
        I = ξᵢ.𝐺
        for (j,ξⱼ) in enumerate(𝓖)
            J = ξⱼ.𝐺
            for ii in 1:3
                for jj in 1:3
                    k[I,J] += gᵢ[i,ii]*G⁻¹[ii,jj]*gⱼ[j,jj]
                end
            end
        end
    end
end

rk = rank(k)