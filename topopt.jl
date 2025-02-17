function design_variables_init()

    # Define the unfiltered design design_variables

    p_reffe = ReferenceFE(lagrangian, Float64, 0)
    Q = TestFESpace(Ω_d, p_reffe, vector_type = Vector{Float64})
    P = Q
    np = num_free_dofs(P)

    # Define the filtered design design_variables

    pf_reffe = ReferenceFE(lagrangian, Float64, 1)
    Qf = TestFESpace(Ω_d, pf_reffe, vector_type = Vector{Float64})
    Pf = Qf

    # We pack everything in a dictionary

    design_params = (; Q, P, Qf, Pf, np)

    return design_params

end

function Filter(r_f, dΩ_d)

    a_f(r_f, u, v) = r_f^2 * (∇(v) ⋅ ∇(u))

    ξ_init = FEFunction(design_params.P, p0)
    op = AffineFEOperator(design_params.Pf, design_params.Qf) do u, v
        ∫(a_f(r, u, v))dΩ_d + ∫(v * u)dΩ_d, ∫(v * ξ_init)dΩ_d
      end
    ξ_f = solve(op)
    return get_free_dof_values(ξ_f)

end

function Threshold(ξ_f; β, η)
    return ((tanh(β * η) + tanh(β * (ξ_f - η))) / (tanh(β * η) + tanh(β * (1.0 - η))))
end

function mat_interpolation_perm(ξ, ε₀, ε₁, type="linear")
    if type == "linear"
        return ε₀ + ξ*(ε1 - ε0)
    elseif type == "quadratic"
        error("Unknown interpolation type")
    end  
end



