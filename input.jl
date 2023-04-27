function import_quad(filename::String)
    elms,nds = ApproxOperator.importmsh(filename)
    nₚ = length(nds)
    nodes = Node{(:𝐼,),1}[]
    data = Dict([:x=>(1,zeros(nₚ)),:y=>(1,zeros(nₚ)),:z=>(1,zeros(nₚ))])
    for (i,p) in enumerate(nds)
        node = Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
    end

    elements = Dict(["Ω"=>Element{:Quad}[],"Ωᵛ"=>Element{:Quad}[],"Γᵍ"=>Element{:Seg2}[],"Γᵗ"=>Element{:Seg2}[]])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 4
    gauss_scheme = :QuadGI4
    nₑ = length(elms["Ω"])

    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*4)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*4)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*4)),
    ])
    for (C,a) in enumerate(elms["Ω"])
        element = Element{:Quad}((c,4,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 4

        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = ApproxOperator.get𝐽(a,ξ,η)*x.w
            push!(𝓖,x)
            s += 4
        end
        g += ng
        push!(elements["Ω"],element)
    end
    
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 1
    gauss_scheme = :QuadGI1
    nₑ = length(elms["Ω"])

    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*4)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*4)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*4)),
    ])
    for (C,a) in enumerate(elms["Ω"])
        element = Element{:Quad}((c,4,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 4

        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = ApproxOperator.get𝐽(a,ξ,η)*x.w
            push!(𝓖,x)
            s += 4
        end
        g += ng
        push!(elements["Ωᵛ"],element)
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 2 
    gauss_scheme = :SegGI2
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵍ"])

    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*2)),
    ])
    for (C,a) in enumerate(elms["Γᵍ"])
        element = Element{:Seg2}((c,2,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 2
       
        𝐿 = ApproxOperator.get𝐿(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w
            push!(𝓖,x)
            s += 2
        end
        g += ng
        push!(elements["Γᵍ"],element)
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 2 
    gauss_scheme = :SegGI2
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵗ"])

    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :n₁=>(2,zeros(ng*nₑ)),
        :n₂=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*2)),
    ])
    for (C,a) in enumerate(elms["Γᵗ"])
        element = Element{:Seg2}((c,2,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 2
       
        𝐿 = ApproxOperator.get𝐿(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w
            push!(𝓖,x)
            s += 2
        end
        g += ng
        push!(elements["Γᵗ"],element)
    end
    return elements,nodes
end

function import_gauss_quadratic(filename::String,s::Symbol)
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
    scheme = ApproxOperator.quadraturerule(s)
    n𝒑 = 21

    elements = Dict([
        "Ω"=>ReproducingKernel{parameters...,:Tri3}[],
        "Γᵗ"=>ReproducingKernel{parameters...,:Seg2}[],
        "Γᵍ"=>ReproducingKernel{parameters...,:Seg2}[]
    ])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = length(scheme[:w])
    ns = 0
    nₑ = length(elms["Ω"])

    for (C,a) in enumerate(elms["Ω"])
        indices = Set{Int}()
        for i in 1:ng
            ξ = scheme[:ξ][i]
            η = scheme[:η][i]
            x,y,z = a(ξ,η)
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Tri3}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Ω"],element)

        c += nc
        g += ng
        ns += nc*ng
    end

    data = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(g)),
        :y=>(2,zeros(g)),
        :z=>(2,zeros(g)),
        :𝑤=>(2,zeros(g)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
    ])
    
    G = 0
    s = 0
    for (C,a) in enumerate(elms["Ω"])
        𝐴 = ApproxOperator.get𝐴(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖,x)
            s += getfield(elements["Ω"][C],:𝓒)[2]
        end
    end
    
    if haskey(elms,"Γᵗ")
        𝓒 = Node{(:𝐼,),1}[]
        𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
        c = 0
        g = 0
        ng = 3
        ns = 0
        gauss_scheme = :SegGI3
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
    gauss_scheme = :SegGI3
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
    
function import_rkgsi_mix(filename1::String,filename2::String)
    elms,nds = ApproxOperator.importmsh(filename1)
    ~,pis = ApproxOperator.importmsh(filename2)
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
    nₚ_𝑝 = length(nds)
    nodes_𝑝 = Node{(:𝐼,),1}[]
    x_𝑝 = zeros(nₚ_𝑝)
    y_𝑝 = zeros(nₚ_𝑝)
    z_𝑝 = zeros(nₚ_𝑝)
    data_𝑝 = Dict([:x=>(1,x_𝑝),:y=>(1,y_𝑝),:z=>(1,z_𝑝)])
    for (i,p) in enumerate(pis)
        node = Node{(:𝐼,),1}((i,),data_𝑝)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes_𝑝,node)
    end
    sp = ApproxOperator.RegularGrid(x,y,z,n=3,γ=5)
    sp_𝑝 = ApproxOperator.RegularGrid(x_𝑝,y_𝑝,z_𝑝,n=3,γ=5)

    parameters = (:Quadratic2D,:□,:CubicSpline)
    scheme_Ω = ApproxOperator.quadraturerule(:TriRK6)
    scheme_Ω̃ = ApproxOperator.quadraturerule(:TriGI3)
    n𝒑 = 21
    n𝒑̃ = 6

    elements = Dict([
        "Ω"=>ReproducingKernel{parameters...,:Tri3}[],
        "Ωᵖ"=>ReproducingKernel{parameters...,:Tri3}[],
        "Ω̃"=>RKGradientSmoothing{parameters...,:Tri3}[],
        "Ω̃ᵖ"=>GRKGradientSmoothing{parameters...,:Tri3}[],
        "Γᵗ"=>ReproducingKernel{parameters...,:Seg2}[],
        "Γᵍ"=>ReproducingKernel{parameters...,:Seg2}[]
    ])

    𝓒 = Node{(:𝐼,),1}[]
    𝓒ᵖ = Node{(:𝐼,),1}[]
    𝓖_Ω = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    𝓖_Ω̃ = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    𝓖_Ωᵖ = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    𝓖_Ω̃ᵖ = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    cᵖ = 0
    g_Ω = 0
    g_Ω̃ = 0
    ng_Ω = 6
    ng_Ω̃ = 3
    ns_Ω = 0
    ns_Ω̃ = 0
    ns_Ωᵖ = 0
    ns_Ω̃ᵖ = 0
    nₑ = length(elms["Ω"])

    𝗚 = zeros(nₚ_𝑝,nₚ_𝑝)
    𝗴₁ = zeros(nₚ_𝑝,nₚ)
    𝗴₂ = zeros(nₚ_𝑝,nₚ)
    for (C,a) in enumerate(elms["Ω"])
        indices = Set{Int}()
        indices_𝑝 = Set{Int}()
        for i in 1:ng_Ω
            ξ = scheme_Ω[:ξ][i]
            η = scheme_Ω[:η][i]
            x,y,z = a(ξ,η)
            union!(indices,sp(x,y,z))
            union!(indices_𝑝,sp_𝑝(x,y,z))
        end
        nc = length(indices)
        nc_Ω̄ = length(indices_𝑝)
        for i in indices
            push!(𝓒,nodes[i])
        end
        for i in indices_𝑝
            push!(𝓒_Ω̄,nodes[i])
        end
        element_Ω = ReproducingKernel{parameters...,:Tri3}((c,nc,𝓒),(g_Ω,ng_Ω,𝓖_Ω))
        element_Ω̃ = RKGradientSmoothing{parameters...,:Tri3}((c,nc,𝓒),(g_Ω̃,ng_Ω̃,𝓖_Ω̃),(g_Ω,ng_Ω,𝓖_Ω))
        element_Ω̄ = GRKGradientSmoothing{parameters...,:Tri3}((c,nc,𝓒),(c_Ω̄,nc_Ω̄,𝓒_Ω̄),(g_Ω̃,ng_Ω̃,𝓖_Ω̄),(g_Ω,ng_Ω,𝓖_Ω),(g_Ω,ng_Ω,𝓖_Ω̂),𝗚,𝗴₁,𝗴₂)
        push!(elements["Ω"],element_Ω)
        push!(elements["Ω̃"],element_Ω̃)
        push!(elements["Ω̄"],element_Ω̄)

        c += nc
        c_Ω̄ += nc_Ω̄
        g_Ω += ng_Ω
        g_Ω̃ += ng_Ω̃
        ns_Ω += nc*ng_Ω
        ns_Ω̃ += nc*ng_Ω̃
        ns_Ω̄ += nc_Ω̄*ng_Ω̃
        ns_Ω̂ += nc_Ω̄*ng_Ω
    end

    D₁ = zeros(g_Ω)
    D₂ = zeros(g_Ω)
    x = zeros(g_Ω)
    y = zeros(g_Ω)
    z = zeros(g_Ω)
    𝑤 = zeros(g_Ω)
    x̃ = zeros(g_Ω̃)
    ỹ = zeros(g_Ω̃)
    z̃ = zeros(g_Ω̃)
    𝑤̃ = zeros(g_Ω̃)
    𝐴 = zeros(nₑ)
    D₁₁ = zeros(nₑ)
    D₁₂ = zeros(nₑ)
    D₂₁ = zeros(nₑ)
    D₂₂ = zeros(nₑ)
    D₃₁ = zeros(nₑ)
    D₃₂ = zeros(nₑ)
    data_𝓖_Ω = Dict([
        :ξ=>(1,scheme_Ω[:ξ]),
        :η=>(1,scheme_Ω[:η]),
        :w=>(1,scheme_Ω[:w]),
        :wᵇ=>(1,scheme_Ω[:wᵇ]),
        :D₁=>(2,D₁),
        :D₂=>(2,D₂),
        :x=>(2,x),
        :y=>(2,y),
        :z=>(2,z),
        :𝑤=>(2,𝑤),
        :𝝭=>(4,zeros(ns_Ω)),
        :𝗠=>(0,zeros(n𝒑)),
    ])
    data_𝓖_Ω̃ = Dict([
        :ξ=>(1,scheme_Ω̃[:ξ]),
        :η=>(1,scheme_Ω̃[:η]),
        :w=>(1,scheme_Ω̃[:w]),
        :x=>(2,x̃),
        :y=>(2,ỹ),
        :z=>(2,z̃),
        :𝑤=>(2,𝑤̃),
        :𝐴=>(3,𝐴),
        :D₁₁=>(3,D₁₁),
        :D₁₂=>(3,D₁₂),
        :D₂₁=>(3,D₂₁),
        :D₂₂=>(3,D₂₂),
        :D₃₁=>(3,D₃₁),
        :D₃₂=>(3,D₃₂),
        :∂𝝭∂x=>(4,zeros(ns_Ω̃)),
        :∂𝝭∂y=>(4,zeros(ns_Ω̃)),
        :∇̃=>(0,zeros(n𝒑̃)),
    ])
    data_𝓖_Ω̄ = Dict([
        :ξ=>(1,scheme_Ω̃[:ξ]),
        :η=>(1,scheme_Ω̃[:η]),
        :w=>(1,scheme_Ω̃[:w]),
        :x=>(2,x̃),
        :y=>(2,ỹ),
        :z=>(2,z̃),
        :𝑤=>(2,𝑤̃),
        :𝐴=>(3,𝐴),
        :D₁₁=>(3,D₁₁),
        :D₁₂=>(3,D₁₂),
        :D₂₁=>(3,D₂₁),
        :D₂₂=>(3,D₂₂),
        :D₃₁=>(3,D₃₁),
        :D₃₂=>(3,D₃₂),
        :𝝭=>(4,zeros(ns_Ω̄)),
        :∂𝝭∂x=>(4,zeros(ns_Ω̄)),
        :∂𝝭∂y=>(4,zeros(ns_Ω̄)),
        :∇̃=>(0,zeros(n𝒑̃)),
    ])
    data_𝓖_Ω̂ = Dict([
        :ξ=>(1,scheme_Ω[:ξ]),
        :η=>(1,scheme_Ω[:η]),
        :w=>(1,scheme_Ω[:w]),
        :wᵇ=>(1,scheme_Ω[:wᵇ]),
        :D₁=>(2,D₁),
        :D₂=>(2,D₂),
        :x=>(2,x),
        :y=>(2,y),
        :z=>(2,z),
        :𝑤=>(2,𝑤),
        :𝝭=>(4,zeros(ns_Ω̂)),
        :∂𝝭∂x=>(4,zeros(ns_Ω̂)),
        :∂𝝭∂y=>(4,zeros(ns_Ω̂)),
        :𝗠=>(0,zeros(n𝒑)),
    ])
    
    G_Ω = 0
    s_Ω = 0
    s_Ω̂ = 0
    s_Ω̄ = 0
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
            x_𝑝 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G_Ω,C,s_Ω̂),data_𝓖_Ω̂)
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
            push!(𝓖_Ω̂,x_𝑝)
            s_Ω += getfield(elements["Ω"][C],:𝓒)[2]
            s_Ω̂ += getfield(elements["Ω̄"][C],:𝓒ᵖ)[2]
        end
        for i in 1:ng_Ω̃
            G_Ω̃ += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G_Ω̃,C,s_Ω̃),data_𝓖_Ω̃)
            x_𝑝 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G_Ω̃,C,s_Ω̄),data_𝓖_Ω̄)
            ξ = x.ξ
            η = x.η
                
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖_Ω̃,x)
            push!(𝓖_Ω̄,x_𝑝)
            s_Ω̃ += getfield(elements["Ω"][C],:𝓒)[2]
            s_Ω̄ += getfield(elements["Ω̄"][C],:𝓒ᵖ)[2]
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
    
