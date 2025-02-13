using Gridap
using Gridap.Geometry

function fea_init(model, order, degree, neumann_tags, source_tags)

    # Define the FEM space

    reffe = ReferenceFE(lagrangian,Float64,order)
    V = TestFESpace(model,reffe, dirichlet_tags=neumann_tags, vector_type=Vector{ComplexF64})
    U = V # mathematically equivalent to TrialFESpace(V,0)
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    Γ_n = BoundaryTriangulation(model;tags=neumann_tags)
    dΓ_n = Measure(Γ_n,degree)

    if source_tags == "Source"
        Γ_s = BoundaryTriangulation(model;tags=source_tags)
        dΓ_s = Measure(Γ_s,degree)

        return U, V, Ω, dΩ, Γ_n, dΓ_n, Γ_s, dΓ_s
    else

        return U, V, Ω, dΩ, Γ_n, dΓ_n, false, false
    end
    
end

function set_tags(model, ε₁, ε₀)

    labels = get_face_labeling(model)
    dimension = num_cell_dims(model)
    tags = get_face_tag(labels,dimension)
    #design_tag = get_tag_from_name(labels,"Design")
    passive_tag = get_tag_from_name(labels,"Passive")

    function ξ(tag)
        #if tag == design_tag
        #   return 0.5*ε₁
        if tag == passive_tag
            return ε₁
        else
            return ε₀
        end
    end

    τ = CellField(tags,Ω) 

    return ξ, τ

end