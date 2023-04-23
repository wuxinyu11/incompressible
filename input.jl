
function import_rkgsi(filename::String)
    elms,nds = ApproxOperator.importmsh(filename)
    nₚ = length(nds)
    nodes = Node{(:𝐼,),1}[]
    x = zeros(nₚ)
    y = zeros(nₚ)
    z = zeros(nₚ)
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    for (i,p) in enumerate(nds)
        node = Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
    end
    sp = ApproxOperator.RegularGrid(x,y,z,n=3,γ=5)

    parameters = (:Quadratic2D,:□,:CubicSpline)
    scheme_Ω = ApproxOperator.quadraturerule(:TriRK6)
    scheme_Ω̃ = ApproxOperator.quadraturerule(:TriGI3)
    n𝒑 = 21
    n𝒑̃ = 6

    elements = Dict([
        "Ω"=>ReproducingKernel{parameters...,:Tri3}[],
        "Ω̃"=>RKGradientSmoothing{parameters...,:Tri3}[],
        "Γᵗ"=>ReproducingKernel{parameters...,:Seg2}[],
        "Γᵍ"=>ReproducingKernel{parameters...,:Seg2}[]
    ])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖_Ω = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    𝓖_Ω̃ = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g_Ω = 0
    g_Ω̃ = 0
    ng_Ω = 6
    ng_Ω̃ = 3
    ns_Ω = 0
    ns_Ω̃ = 0
    nₑ = length(elms["Ω"])

    for (C,a) in enumerate(elms["Ω"])
        indices = Set{Int}()
        for i in 1:ng_Ω
            ξ = scheme_Ω[:ξ][i]
            η = scheme_Ω[:η][i]
            x,y,z = a(ξ,η)
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element_Ω = ReproducingKernel{parameters...,:Tri3}((c,nc,𝓒),(g_Ω,ng_Ω,𝓖_Ω))
        element_Ω̃ = RKGradientSmoothing{parameters...,:Tri3}((c,nc,𝓒),(g_Ω̃,ng_Ω̃,𝓖_Ω̃),(g_Ω,ng_Ω,𝓖_Ω))
        push!(elements["Ω"],element_Ω)
        push!(elements["Ω̃"],element_Ω̃)

        c += nc
        g_Ω += ng_Ω
        g_Ω̃ += ng_Ω̃
        ns_Ω += nc*ng_Ω
        ns_Ω̃ += nc*ng_Ω̃
    end

    data_𝓖_Ω = Dict([
        :ξ=>(1,scheme_Ω[:ξ]),
        :η=>(1,scheme_Ω[:η]),
        :w=>(1,scheme_Ω[:w]),
        :wᵇ=>(1,scheme_Ω[:wᵇ]),
        :D₁=>(2,zeros(g_Ω)),
        :D₂=>(2,zeros(g_Ω)),
        :x=>(2,zeros(g_Ω)),
        :y=>(2,zeros(g_Ω)),
        :z=>(2,zeros(g_Ω)),
        :𝑤=>(2,zeros(g_Ω)),
        :𝝭=>(4,zeros(ns_Ω)),
        :𝗠=>(0,zeros(n𝒑)),
    ])
    data_𝓖_Ω̃ = Dict([
        :ξ=>(1,scheme_Ω̃[:ξ]),
        :η=>(1,scheme_Ω̃[:η]),
        :w=>(1,scheme_Ω̃[:w]),
        :x=>(2,zeros(g_Ω̃)),
        :y=>(2,zeros(g_Ω̃)),
        :z=>(2,zeros(g_Ω̃)),
        :𝑤=>(2,zeros(g_Ω̃)),
        :𝐴=>(3,zeros(nₑ)),
        :D₁₁=>(3,zeros(nₑ)),
        :D₁₂=>(3,zeros(nₑ)),
        :D₂₁=>(3,zeros(nₑ)),
        :D₂₂=>(3,zeros(nₑ)),
        :D₃₁=>(3,zeros(nₑ)),
        :D₃₂=>(3,zeros(nₑ)),
        :∂𝝭∂x=>(4,zeros(ns_Ω̃)),
        :∂𝝭∂y=>(4,zeros(ns_Ω̃)),
        :∇̃=>(0,zeros(n𝒑̃)),
    ])
    
    G_Ω = 0
    s_Ω = 0
    G_Ω̃ = 0
    s_Ω̃ = 0
    for (C,a) in enumerate(elms["Ω"])
        𝐴 = ApproxOperator.get𝐴(a)
        x₁ = a.vertices[1].x
        x₂ = a.vertices[2].x
        x₃ = a.vertices[3].x
        y₁ = a.vertices[1].y
        y₂ = a.vertices[2].y
        y₃ = a.vertices[3].y
        D₁₁ = y₃-y₂
        D₁₂ = x₂-x₃
        D₂₁ = y₁-y₃
        D₂₂ = x₃-x₁
        D₃₁ = y₂-y₁
        D₃₂ = x₁-x₂
        for i in 1:ng_Ω
            G_Ω += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G_Ω,C,s_Ω),data_𝓖_Ω)
            ξ = x.ξ
            η = x.η
            D₁ = 0.0
            D₂ = 0.0
            if ξ ≈ 0.0 (D₁ += D₁₁;D₂ += D₁₂) end
            if η ≈ 0.0 (D₁ += D₂₁;D₂ += D₂₂) end
            if ξ+η ≈ 1.0 (D₁ += D₃₁;D₂ += D₃₂) end
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            x.D₁ = D₁
            x.D₂ = D₂
            push!(𝓖_Ω,x)
            s_Ω += getfield(elements["Ω"][C],:𝓒)[2]
        end
        for i in 1:ng_Ω̃
            G_Ω̃ += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G_Ω̃,C,s_Ω̃),data_𝓖_Ω̃)
            ξ = x.ξ
            η = x.η
                
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖_Ω̃,x)
            s_Ω̃ += getfield(elements["Ω"][C],:𝓒)[2]
        end
        elements["Ω̃"][C].𝐴 = 𝐴
        elements["Ω̃"][C].D₁₁ = D₁₁
        elements["Ω̃"][C].D₁₂ = D₁₂
        elements["Ω̃"][C].D₂₁ = D₂₁
        elements["Ω̃"][C].D₂₂ = D₂₂
        elements["Ω̃"][C].D₃₁ = D₃₁
        elements["Ω̃"][C].D₃₂ = D₃₂
    end
    
    if haskey(elms,"Γᵗ")
        𝓒 = Node{(:𝐼,),1}[]
        𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
        c = 0
        g = 0
        ng = 3
        ns = 0
        gauss_scheme = :SegRK3
        scheme = ApproxOperator.quadraturerule(gauss_scheme)
        nₑ = length(elms["Γᵗ"])

        for (C,a) in enumerate(elms["Γᵗ"])
            indices = Set{Int}()
            for i in 1:ng
                ξ = scheme[:ξ][i]
                x,y,z = a(ξ)
                union!(indices,sp(x,y,z))
            end
            nc = length(indices)
            for i in indices
                push!(𝓒,nodes[i])
            end
            element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
            push!(elements["Γᵗ"],element)
            c += nc
            g += ng
            ns += ng*nc
        end

        G = 0
        s = 0
        data_𝓖 = Dict([
            :ξ=>(1,scheme[:ξ]),
            :w=>(1,scheme[:w]),
            :x=>(2,zeros(ng*nₑ)),
            :y=>(2,zeros(ng*nₑ)),
            :z=>(2,zeros(ng*nₑ)),
            :𝑤=>(2,zeros(ng*nₑ)),
            :n₁=>(3,zeros(nₑ)),
            :n₂=>(3,zeros(nₑ)),
            :𝗠=>(0,zeros(n𝒑)),
            :𝝭=>(4,zeros(ns))
        ])
        for (C,a) in enumerate(elms["Γᵗ"])
            𝐿 = ApproxOperator.get𝐿(a)
            x₁ = a.vertices[1].x
            x₂ = a.vertices[2].x
            y₁ = a.vertices[1].y
            y₂ = a.vertices[2].y
            n₁ = (y₂-y₁)/𝐿
            n₂ = (x₁-x₂)/𝐿
            for i in 1:ng
                G += 1
                x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
                ξ = x.ξ
                x_,y_,z_ = a(ξ)
                x.x = x_
                x.y = y_
                x.z = z_
                x.𝑤 = 𝐿*x.w/2
                push!(𝓖,x)
                s += getfield(elements["Γᵗ"][C],:𝓒)[2]
            end
            elements["Γᵗ"][C].n₁ = n₁
            elements["Γᵗ"][C].n₂ = n₂
        end
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = 3
    ns = 0
    gauss_scheme = :SegRK3
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵍ"])

    for (C,a) in enumerate(elms["Γᵍ"])
        indices = Set{Int}()
        for i in 1:ng
            ξ = scheme[:ξ][i]
            x,y,z = a(ξ)
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Γᵍ"],element)
        c += nc
        g += ng
        ns += ng*nc
    end
       
    G = 0
    s = 0
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :n₁=>(3,zeros(nₑ)),
        :n₂=>(3,zeros(nₑ)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns))
    ])
    for (C,a) in enumerate(elms["Γᵍ"])
        𝐿 = ApproxOperator.get𝐿(a)
        x₁ = a.vertices[1].x
        x₂ = a.vertices[2].x
        y₁ = a.vertices[1].y
        y₂ = a.vertices[2].y
        n₁ = (y₂-y₁)/𝐿
        n₂ = (x₁-x₂)/𝐿
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w/2
            push!(𝓖,x)
            s += getfield(elements["Γᵍ"][C],:𝓒)[2]
        end
        elements["Γᵍ"][C].n₁ = n₁
        elements["Γᵍ"][C].n₂ = n₂
    end


    return elements,nodes
end
    
