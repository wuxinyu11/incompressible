using Revise,ApproxOperator,TOML,CairoMakie

file_nodes = "./msh/square_node_12.msh"
file_elements = "./msh/square_bubble_24.msh"
config = TOML.parsefile("./toml/plot_msh.toml")
elements,nodes = importmsh(file_elements,file_nodes,config)

f = Figure()

# axis
ax = Axis(f[1, 1])
limits!(ax, -0.1,1.1,-0.1,1.1) 
# limits!(ax, 0,1,0,1) 
hidespines!(ax)
hidedecorations!(ax)
ax.aspect = AxisAspect(1)


# nodes
data = getfield(nodes[1],:data)
x = data[:x][2]
y = data[:y][2]
scatter!(x,y, 
    marker=:circle,
    markersize = 20,
    color = :black
)

# boundaries
elements["Γ"] = elements["Γᵗ"]∪elements["Γᵍ"]
for elm in elements["Γ"]
    x = [x.x for x in elm.𝓒]
    y = [x.y for x in elm.𝓒]
    lines!(x,y,linewidth = 3, color = :black)
end

# elements
for elm in elements["Ω"]
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]
    lines!(x,y,linestyle = :dash, linewidth = 1.5, color = :black)
end

f