using Gridap.Fields

function s_PML(x,σ,k,LH,d_pml)
    u = abs.(Tuple(x)).-LH./2  # get the depth into PML
    return @. ifelse(u > 0,  1+(1im*σ/k)*(u/d_pml)^2, $(1.0+0im))
end

function ds_PML(x,σ,k,LH,d_pml)
    u = abs.(Tuple(x)).-LH./2 # get the depth into PML
    ds = @. ifelse(u > 0, (2im*σ/k)*(1/d_pml)^2*u, $(0.0+0im))
    return ds.*sign.(Tuple(x))
end

struct Λ<:Function
    σ::Float64
    k::Float64
    LH::NTuple{2,Float64}
    d_pml::Float64
end

function (Λf::Λ)(x)
    s_x,s_y = s_PML(x,Λf.σ,Λf.k,Λf.LH,Λf.d_pml)
    return VectorValue(1/s_x,1/s_y)
end

Fields.∇(Λf::Λ) = x->TensorValue{2,2,ComplexF64}(-(Λf(x)[1])^2*ds_PML(x,Λf.σ,Λf.k,Λf.LH,Λf.d_pml)[1],0,0,-(Λf(x)[2])^2*ds_PML(x,Λf.σ,Λf.k,Λf.LH,Λf.d_pml)[2])

function pml(L,H,k,d_pml)

    Rpml = 1e-12      # Tolerance for PML reflection
    σ = -3/4*log(Rpml)/d_pml # σ_0
    LH = (L,H)
    Λf = Λ(σ,k,LH,d_pml)

    return Λf

end


