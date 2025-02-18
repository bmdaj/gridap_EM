using ChainRulesCore, Zygote
import ChainRulesCore: rrule
NO_FIELDS = ZeroTangent()

dξt_dξf(ξ_f, β, η) = β * (1.0 - tanh(β * (ξ_f - η))^2) / (tanh(β * η) + tanh(β * (1.0 - η)))
dε_dξf(ξ_f, ε₀, ε₁, β, η) = (ε₁ - ε₀) * dξt_dξf(ξ_f, β, η)
dA_dξf(u, v, ξ_th; β, η) = -(k^2 * (((ξ -> dε_dξf(ξ, ε₀, ε₁, β, η)) ∘ ξ_th) * v * u))

function ξf_ξ(ξ; r_f, dΩ_d)
    pf_vec = Filter(ξ; r_f, dΩ_d)
    pf_vec
end

function rrule(::typeof(ξf_ξ), ξ; r_f, dΩ_d)
    function ξ_f_pullback(dfom_dξf)
      NO_FIELDS, dξf_dξ(dfom_dξf; r_f, dΩ_d, design_params)
    end
    ξf_ξ(ξ; r_f, dΩ_d), ξ_f_pullback
  end

function dξf_dξ(dfom_dξf; r_f, dΩ_d, design_params)
    Af = assemble_matrix(design_params.Pf, design_params.Qf) do u, v
        ∫(a_f(r_f, u, v))dΩ_d + ∫(v * u)dΩ_d
    end
    λ_adj_vec = Af' \ dfom_dξf
    λ_adjh= FEFunction(design_params.Pf, λ_adj_vec)
    l_temp(dξ) = ∫(λ_adjh * dξ)dΩ_d
    return assemble_vector(l_temp, design_params.P)
end

function FOM_ξf(ξ::Vector; U, V, r_f, β, η, dΩ, dΩ_d, dΓ_s, design_params)
    ξ_f_vec = ξf_ξ(ξ; r_f, dΩ_d)
    FOM_eval(ξ_f_vec; U, V, β, η, dΩ, dΓ_s, design_params)
end

function dFOM_dξ(ξ::Vector, grad::Vector; U, V, r_f, β, η, dΩ, dΩ_d, dΓ_s, design_params)
    if length(grad) > 0
        dfom_dξ, = Zygote.gradient(x -> FOM_ξf(x; U, V, r_f, β, η, dΩ, dΩ_d, dΓ_s, design_params), ξ)
        grad[:] = dfom_dξ
    end
    FOM_value = FOM_ξf(ξ::Vector; U, V, r_f, β, η, dΩ, dΩ_d, dΓ_s, design_params)
    open("gvalue.txt", "a") do io
        write(io, "$FOM_value \n")
    end
    FOM_value
end