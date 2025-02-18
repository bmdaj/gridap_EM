function design_variables_init()

    # Define the unfiltered design design_variables

    ξ_reffe = ReferenceFE(lagrangian, Float64, 0)
    Q = TestFESpace(Ω_d, ξ_reffe, vector_type = Vector{Float64})
    P = Q
    np = num_free_dofs(P)

    # Define the filtered design design_variables

    ξf_reffe = ReferenceFE(lagrangian, Float64, 1)
    Qf = TestFESpace(Ω_d, ξf_reffe, vector_type = Vector{Float64})
    Pf = Qf

    # We pack everything in a dictionary

    design_params = (; Q, P, Qf, Pf, np)

    return design_params

end

a_f(r_f, u, v) = r_f^2 * (∇(v) ⋅ ∇(u))

function Filter(ξ; r_f, dΩ_d)

    ξ_init = FEFunction(design_params.P, ξ)
    op = AffineFEOperator(design_params.Pf, design_params.Qf) do u, v
        ∫(a_f(r_f, u, v))dΩ_d + ∫(v * u)dΩ_d, ∫(v * ξ_init)dΩ_d
      end
    ξ_f = solve(op)

    return get_free_dof_values(ξ_f)

end

function Threshold(ξ_f; β, η)
    return ((tanh(β * η) + tanh(β * (ξ_f - η))) / (tanh(β * η) + tanh(β * (1.0 - η))))
end

function ξ_inter(ξ_th, ε₀, ε₁)
    return ε₀ + ξ_th*(ε₁ - ε₀) 
end



