
using  YAML, ApproxOperator, XLSX, LinearAlgebra 

r = 4

ndiv2 = 96
index = [12,24,48,96]

ops = [
    Operator{:∫∇v∇udΩ}(:k=>1.0),
    Operator{:∫vbdΩ}(),
    Operator{:∫vtdΓ}(),
    Operator{:∫∇𝑛vgdΓ}(:k=>1.0),
    Operator{:∫∇̄𝑛vgdΓ}(:k=>1.0),
    Operator{:H₁}()
]

 
for n in 22:96
ndiv1 = Int(ndiv2*n/96)
# n =24
#    ndiv1 = 24

path1 = "./msh/heat_"*string(ndiv1)*".msh"
path2 = "./msh/heat_"*string(ndiv2)*".msh"

# config = YAML.load_file("./yml/heat_quadratic.yml")
config = YAML.load_file("./yml/heat_cubic.yml")

elements, nodes = importmsh(path1,path2,config)

set_memory_𝗠!(elements["Ω̃"],:∇̃)
set_memory_𝗠!(elements["Γᵍ"],:𝝭,:∇̃)
nₚ = length(nodes)
# s = 2.5/ndiv2*ones(nₚ)
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
L2=log10(l2)

# l2 = 0.
# h1 = 0.
# if det(k) ≈ 0.
#     l2 = NaN
# else
#     d = k\f

#     push!(nodes,:d=>d)
#     set∇𝝭!(elements["Ωᴳ"])
#     prescribe!(elements["Ωᴳ"],:u=>(x,y,z)->(x+y)^r)
#     prescribe!(elements["Ωᴳ"],:∂u∂x=>(x,y,z)->r*(x+y)^abs(r-1))
#     prescribe!(elements["Ωᴳ"],:∂u∂y=>(x,y,z)->r*(x+y)^abs(r-1))
#     h1,l2 = ops[6](elements["Ωᴳ"])
# end


#     push!(nodes,:d=>d)
#     set∇𝝭!(elements["Ωᴳ"])
#     prescribe!(elements["Ωᴳ"],:u=>(x,y,z)->(x+y)^r)
#     h1,l2 = ops[6](elements["Ωᴳ"])
# end

# h = log10(1/ndiv1)

# f = check∇𝝭(elements["Ω"])
# f = check∇𝝭(elements["Ω̃"])
# f = check∇𝝭(elements["Γᵍ"])
XLSX.openxlsx("./xlsx/heat.xlsx", mode="rw") do xf
    # row = Char(64+findfirst(n_->n_==n,1:ndiv2))
    row = Char(64+findfirst(n_->n_==ndiv2,index))
    𝐿₂ = xf[2]
    𝐻₁ = xf[3]
    # ind = findfirst(n_->n_==ndiv2,index)+1
    ind = findfirst(n_->n_==n,1:ndiv2)
    row = row*string(ind)
    𝐿₂[row] = log10(l2)
    𝐻₁[row] = log10(h1)
end
end