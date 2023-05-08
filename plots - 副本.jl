using YAML,ApproxOperator,Plots, LinearAlgebra

## 1D style
# linestyle = (:solid,2)
# markstyle = (:circle,6,Plots.stroke(0, :black))
# xlims_ = (0,1)
# ylims_ = (-1,1)

## 2D style
# patch_test:
#   linewidth: 2->1.5, 12->1
#   markersize: 2->6, 6->5

# xlims_ = (-0.05,1.05)
# ylims_ = (-0.05,1.05)

# plate_with_hole:
#   linewidth: 3,6->2, 12->1, 24->0.5
#   markersize: 3->6, 6->5, 12->3.5, 24->1

# xlims_ = (-0.1,5.1)
# ylims_ = (-0.1,5.1)
# xlims_ = (-0.1,5.1)
# ylims_ = (-0.1,5.1)

# cantilever:
#   linewidth: 4->2, 8->1.5, 16->1, 32->0.5
#   markersize: 4->5, 8->4, 16->3, 32->1
# xlims_ = (-0.5,48.5)
# ylims_ = (-6.5,6.5)
#  cook:
#   linewidth: 4->2, 8->1.5, 16->1,  20->0.6  32->0.5
#   markersize: 4->5, 8->4, 16->3, 20->2 32->1
xlims_ = (0,48.0)
ylims_ = (0,60.0)

linestyle = (:solid,2)
markstyle = (:circle,5,Plots.stroke(0, :black))

function plotmesh(as::Vector{T}) where T<:ApproxOperator.AbstractElement
    p = plot(;framestyle=:none,legend = false,background_color=:transparent,aspect_ratio=:equal)
    for a in as
        plotmesh(a,p)
    end
    return p
end

function plotmesh(a::T,p::Plots.Plot) where T<:ApproxOperator.AbstractElement{:Seg2}
    plot!([a.𝓒[i].x for i in 1:2],zeros(2),color=:black,line=linestyle,marker=markstyle)
end

function plotmesh(a::T,p::Plots.Plot) where T<:ApproxOperator.AbstractElement{:Tri3}
    𝓒 = collect(a.𝓒)
    plot!([𝓒[i].x for i in (1,2,3,1)],[𝓒[i].y for i in (1,2,3,1)],color=:black,line=linestyle,marker=markstyle)
end

# elements, nodes =ApproxOperator.importmsh("./msh/beam_10.msh")
# p = plotmesh(elements["Ω"])
# p;savefig(p,"./figure/bar_mesh.svg")

# elements, nodes = ApproxOperator.importmsh("./msh/patchtest.msh",config)
# p = plotmesh(elements["Ω"])
# savefig(p,"./figure/patchtest.svg")

# config = YAML.load_file("./yml/cantilever_gauss_nitsche_quadratic.yml")
# elements, nodes = ApproxOperator.importmsh("./msh/cantilever_10.msh",config)


include("input.jl")
# elements, nodes = import_tri3("./msh/cook_membrance_20.msh")
elements, nodes = import_tri3("./msh/cantilever_4.msh")

# elements, nodes = ApproxOperator.importmsh("./msh/beam_10.msh")
# include("input.jl")
# elements, nodes = import_gauss_quadratic("./msh/cantilever_10.msh",:TriGI3)
# set𝝭!(elements["Ω"])
p = plotmesh(elements["Ω"])

# lw = 1.5
# plot!([-10/3,20/3,-10/3,-10/3],[-10/3^0.5,0,10/3^0.5,-10/3^0.5],color=:black,line=(:solid,lw))
# plot!([0,0],[1,2],color=:black,line=(:solid,lw))
# plot!([1,2],[0,0],color=:black,line=(:solid,lw))
# curves!([0,0.552284749831,1,1],[1,1,0.552284749831,0],color=:black,line=(:solid,lw))
# curves!([0,2*0.552284749831,2,2],[2,2,2*0.552284749831,0],color=:black,line=(:solid,lw))
# plot!([nodes[i].x for i in (1,2,3,4,1)],[nodes[i].y for i in (1,2,3,4,1)],color=:black,line=(:solid,2.0))
# savefig(p,"./figure/patch_test.svg")
# savefig(p,"./figure/cantilever_10.svg")

# elements, nodes = importmsh("./msh/plate_with_hole_24.msh")
# p = plotmesh(elements["Ω"])
# savefig(p,"./figure/cook_membrance_20.svg")
savefig(p,"./figure/cantilever_4.svg")
plot(p)
