
using  YAML, ApproxOperator, XLSX, LinearAlgebra 

r = 3
ndiv2 = 12

index = [12,24,48,96]

ops = [
    Operator{:∫∇v∇udΩ}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫∇𝑛vgdΓ}(:k=>1.0),
    Operator{:∫vgdΓ}(:α=>1e3),
    Operator{:L₂}()
]

path2 = "./msh/bar_"*string(ndiv2)*".msh"

# for n in 1:12
# n = 6

# ndiv1 = Int(ndiv2*n/12)
ndiv1=1
path1 = "./msh/bar_"*string(ndiv1)*".msh"

config2 = YAML.load_file("./yml/bar1.yml")
elements, nodes = importmsh(path1,path2,config2)

nₚ = length(nodes)

set_memory_𝗠!(elements["Ω̃"],:∇̃)

s = 2.5/ndiv2*ones(nₚ)

push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
set∇𝝭!(elements["Γᵍ"])

prescribe!(elements["Ω"],:b=>(x,y,z)->-r*(r-1)*x^abs(r-2))
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->x^r)
prescribe!([elements["Γᵍ"][1]],:n₁=>(x,y,z)->-1.0)
prescribe!(elements["Γᵍ"][2],:n₁=>(x,y,z)->1.0)



k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω̃"],k)
ops[2](elements["Ω"],f)
# ops[3](elements["Γᵗ"],f)
ops[4](elements["Γᵍ"],k,f)
ops[5](elements["Γᵍ"],k,f)

# l2 = 0.
# # if rank(k) < nₚ
# if det(k) ≈ 0.
#     l2 = NaN
# else
    d = k\f

    push!(nodes,:d=>d)
    set𝝭!(elements["Ωᴳ"])
    prescribe!(elements["Ωᴳ"],:u=>(x,y,z)->x^r)
    l2 = ops[6](elements["Ωᴳ"])
L2=log10(l2)


