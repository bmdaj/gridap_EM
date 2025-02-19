using ChainRulesCore, Zygote
import ChainRulesCore: rrule
NO_FIELDS = ZeroTangent()

dξt_dξf(ξ_f, β, η) = β * (1.0 - tanh(β * (ξ_f - η))^2) / (tanh(β * η) + tanh(β * (1.0 - η)))
dε_dξf(ξ_f, ε₀, ε₁, β, η) = (ε₁ - ε₀) * dξt_dξf(ξ_f, β, η)
dA_dξf(u, v, ξ_th; β, η) = -(k^2 * (((ξ -> dε_dξf(ξ, ε₀, ε₁, β, η)) ∘ ξ_th) * v * u))

function rrule(::typeof(FOM_eval), ξ_f_vec; β, η, fem_params, design_params)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * dFOM_dξf(ξ_f_vec; β, η, fem_params, design_params)
    end
    FOM_eval(ξ_f_vec; β, η, fem_params, design_params), U_pullback
end

function dFOM_dξf(ξ_f_vec; β, η, fem_params, design_params)

    ξ_fh = FEFunction(design_params.Pf, ξ_f_vec)
    ξ_th = (ξ_f -> Threshold(ξ_f; β, η)) ∘ ξ_fh
    A_mat = MatrixA(ξ_th; fem_params.U, fem_params.V)
    b_vec = assemble_vector(v->(∫(v)dΓ_s), fem_params.V)
    Ez_vec = A_mat \ b_vec

    O_mat = MatrixFOM(fem_params)

    Ez = FEFunction(fem_params.U, Ez_vec)
    λ_adj_vec =  A_mat' \ (O_mat * Ez_vec)
    λ_adj_conjh = FEFunction(fem_params.U, conj(λ_adj_vec))

    l_temp(dξ) = ∫(real(-2 * dA_dξf(Ez, λ_adj_conjh, ξ_fh; β, η)) * dξ)fem_params.dΩ_d
    dfom_dξf = assemble_vector(l_temp, design_params.Pf)

    return dfom_dξf
end

function ξf_ξ(ξ; r_f, fem_params, design_params)
    pf_vec = Filter(ξ; r_f, fem_params, design_params)
    pf_vec
end

function rrule(::typeof(ξf_ξ), ξ; r_f, fem_params, design_params)
    function ξ_f_pullback(dfom_dξf)
      NO_FIELDS, dξf_dξ(dfom_dξf; r_f, fem_params, design_params)
    end
    ξf_ξ(ξ; r_f, fem_params, design_params), ξ_f_pullback
  end

function dξf_dξ(dfom_dξf; r_f, fem_params, design_params)
    Af = assemble_matrix(design_params.Pf, design_params.Qf) do u, v
        ∫(a_f(r_f, u, v))fem_params.dΩ_d + ∫(v * u)fem_params.dΩ_d
    end
    λ_adj_vec = Af' \ dfom_dξf
    λ_adjh= FEFunction(design_params.Pf, λ_adj_vec)
    l_temp(dξ) = ∫(λ_adjh * dξ)fem_params.dΩ_d
    return assemble_vector(l_temp, design_params.P)
end

function FOM_ξf(ξ::Vector; r_f, β, η, fem_params, design_params)
    ξ_f_vec = ξf_ξ(ξ; r_f, fem_params, design_params)
    FOM_eval(ξ_f_vec; β, η, fem_params, design_params)
end

function dFOM_dξ(ξ::Vector, grad::Vector; r_f, β, η, fem_params, design_params)
    if length(grad) > 0
        dfom_dξ, = Zygote.gradient(x -> FOM_ξf(x; r_f, β, η, fem_params, design_params), ξ)
        grad[:] = dfom_dξ
    end
    FOM_value = FOM_ξf(ξ::Vector; r_f, β, η, fem_params, design_params)
    open("gvalue.txt", "a") do io
        write(io, "$FOM_value \n")
    end
    FOM_value
end