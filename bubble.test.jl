
using  YAML, ApproxOperator, XLSX, LinearAlgebra 

r = 2
ndiv2 = 12

ops = [
    Operator{:∫∇v∇udΩ}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫∇𝑛vgdΓ}(:k=>1.0),
    Operator{:∫∇̄𝑛vgdΓ}(:k=>1.0),
    Operator{:H₁}()
]

 
# ndiv1 = 4

# path1 = "./msh/heat_"*string(ndiv1)*".msh"
 path1 = "./msh/square_bubble_78-3.msh"

path2 = "./msh/heat_"*string(ndiv2)*".msh"

config = YAML.load_file("./yml/heat_quadratic.yml")
# config = YAML.load_file("./yml/heat_cubic.yml")

elements, nodes = importmsh(path1,path2,config)

set_memory_𝗠!(elements["Ω̃"],:∇̃)
set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∇̃)
nₚ = length(nodes)
s = 3.5/ndiv2*ones(nₚ)

push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set𝝭!(elements["Ω∩Γᵍ"])
set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
set∇̃𝝭!(elements["Γᵍ"],elements["Ω∩Γᵍ"])
set𝝭!(elements["Γᵍ"])
set∇̄𝝭!(elements["Γᵍ"])

# set𝒏!(elements["Γᵍ"])
prescribe!(elements["Ω"],:b=>(x,y,z)->-2*r*(r-1)*(x+y)^abs(r-2))
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->(x+y)^r)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)
ops[4](elements["Γᵍ"],k,f)
ops[5](elements["Γᵍ"],k,f)

 d = k\f

push!(nodes,:d=>d)
set∇𝝭!(elements["Ωᴳ"])
prescribe!(elements["Ωᴳ"],:u=>(x,y,z)->(x+y)^r)
prescribe!(elements["Ωᴳ"],:∂u∂x=>(x,y,z)->r*(x+y)^abs(r-1))
prescribe!(elements["Ωᴳ"],:∂u∂y=>(x,y,z)->r*(x+y)^abs(r-1))
h1,l2 = ops[6](elements["Ωᴳ"])
# L2=log10(l2)
