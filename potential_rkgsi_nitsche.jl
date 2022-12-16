using Revise,ApproxOperator,TOML

file_nodes = "./msh/square_node_12.msh"
file_elements = "./msh/square_bubble_24.msh"
config = TOML.parsefile("./toml/potential_2D_rkgsi_nitsche_cubic.toml")

elements,nodes = importmsh(file_elements,file_nodes,config)
# elements,nodes = importmsh(file_nodes,config)
nₚ = length(nodes)
s = 3.5/12*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

# shape functions
set𝝭!(elements["Ω"])
set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
set𝝭!(elements["Γᵗ"])
set∇₂𝝭!(elements["Γᵍ"])

# prescribing
r = 1
u(x,y,z) = (x+y)^r
∂u∂x(x,y,z) = r*(x+y)^abs(r-1)
∂u∂y(x,y,z) = r*(x+y)^abs(r-1)
∂²u∂x²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
∂²u∂y²(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
t(x,y,z,n₁,n₂) = ∂u∂x(x,y,z)*n₁+∂u∂y(x,y,z)*n₂
b(x,y,z) = -(∂²u∂x²(x,y,z)+∂²u∂y²(x,y,z))

prescribe!(elements["Ω"],:b=>b)
prescribe!(elements["Γᵍ"],:g=>u)

# operator
ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫∇𝑛vgds}(:k=>1.0,:α=>1e3),
    Operator{:H₁}()
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γᵗ"],f)
ops[4](elements["Γᵍ"],k,f)

d = k\f

# rank(k)

config = TOML.parsefile("./toml/error_cubic.toml")
elms,nds = importmsh(file_nodes,config)
push!(nds,:d=>d)
push!(nds,:s₁=>s,:s₂=>s,:s₃=>s)
set∇₂𝝭!(elms["Ω"])
set_memory_𝝭!(elms["Ω"],:∂𝝭∂z)
prescribe!(elms["Ω"],:u=>u)
prescribe!(elms["Ω"],:∂u∂x=>∂u∂x)
prescribe!(elms["Ω"],:∂u∂y=>∂u∂y)
H₁,L₂ = ops[5](elms["Ω"])