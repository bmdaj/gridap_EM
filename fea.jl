using Gridap
using Gridap.Geometry

function fea_init(model, order, degree, dirichlet_tags="None", neumann_tags="None", source_tags="None")

    # Define the FEM space

    reffe = ReferenceFE(lagrangian,Float64,order)

    if dirichlet_tags == "None"
        V = TestFESpace(model,reffe, vector_type=Vector{ComplexF64})
    else
        V = TestFESpace(model,reffe, dirichlet_tags=dirichlet_tags, vector_type=Vector{ComplexF64})
    end

    U = V # mathematically equivalent to TrialFESpace(V,0)
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)

    if neumann_tags == "None"
        
        Γ_n = false
        dΓ_n = false
    
    else 

        Γ_n = BoundaryTriangulation(model;tags=neumann_tags)
        dΓ_n = Measure(Γ_n,degree)

    end

   
    if source_tags == "None"

        Γ_s = false
        dΓ_s = false

    else   

        Γ_s = BoundaryTriangulation(model;tags=source_tags)
        dΓ_s = Measure(Γ_s,degree)

        
    end

    return U, V, Ω, dΩ, Γ_n, dΓ_n, Γ_s, dΓ_s    

end

function set_tags(model, εₛ, tag_list, ε₀)

    labels = get_face_labeling(model)
    dimension = num_cell_dims(model)
    tags = get_face_tag(labels,dimension)
    τ = CellField(tags,Ω) 

    function ξ(tag)

        for i in range(1, length(tag_list))

            tag_name = get_tag_from_name(labels,tag_list[i])

            if tag == tag_name
                return εₛ[i]
            end            
        end
        return ε₀
    end

    return ξ, τ

end