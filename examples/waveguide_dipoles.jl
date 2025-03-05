using GridapGmsh
using Gridap.Fields
using GridapEmbedded

include("../fea.jl")
include("../plotter.jl")
include("../pml.jl")

# We stary by defining the parameters of the problem:

print("Defining model parameters...\n")

λ = 35.0                           # Wavelength (arbitrary unit)
k = 2 * π / λ                      # Wave number 
ε₁ = 4.0                           # Relative electric permittivity for the material
ε₀ = 1.0                           # Relative electric permittivity for the background
εₛ = [0.5 * ε₁, ε₁]                # List of relative permittivities
tag_list = ["Design", "Passive"]   # List of tag names
out = true                         # Output the results to a file

w_tot = 400.0
h_tot = 200.0
d_pml = λ                          # Thickness of the PML
# We import the geometry of the problem from the file `geometry.msh`:

print("Importing mesh...\n")

model = GmshDiscreteModel("waveguide_mesh.msh")

print("Outputing mesh...\n")

if out == true
    writevtk(model,"geom/geometry")
end

# FEA setup:

print("Setting up FEA...\n")

order = 1 
degree = 2
dirichlet_tags = "None"
neumann_tags = "Edges"
source_tags = "None"

U, V, Ω, dΩ, Γ_n, dΓ_n, Γ_s, dΓ_s = fea_init(model, order, degree, dirichlet_tags, neumann_tags, source_tags) 

function dipole_source(pol, δ ,ϵ, x_dip, y_dip, x)

    if pol == "X"

        if abs(x[1] - x_dip - δ) < ϵ && abs(x[2] - y_dip) < ϵ
            return 1.0

        elseif abs(x[1] - x_dip + δ) < ϵ && abs(x[2]- y_dip) < ϵ
            return -1.0

        else
            return 0.0
        end
    elseif pol == "Y"

        if abs(x[1] - x_dip) < ϵ && abs(x[2] - y_dip - δ) < ϵ
            return 1.0

        elseif abs(x[1] - x_dip) < ϵ && abs(x[2] - y_dip + δ) < ϵ
            return -1.0

        else
            return 0.0
        end
    end 
    return 0.0
end

δ = 35.0/50.0 # Separation between the dipoles
ϵ = 35.0/50.0 # Mesh tolerance for dipole separation

N = 2

d_off_dip = λ/2.0

x_dip_list = [75.0/2-d_off_dip, 75.0/2+d_off_dip]
y_dip_list = [0.0, 0.0]

source_functions = []

for pol in ["X", "Y"]
    for i in 1:N
        x_dip_it = x_dip_list[i]
        y_dip_it = y_dip_list[i]
        push!(source_functions, x -> dipole_source(pol, δ, ϵ, x_dip_it, y_dip_it, x))
    end
end

ε_tag, τ = set_tags(model, εₛ, tag_list, ε₀)                                         # We set the permittivity tags for the design and passive regions 

Λf = pml(w_tot, h_tot, k, d_pml)                                                                    # We set the PML function

w(x) = -1.0im * k                                                                                   # absorbing boundary condition coefficient
a(u, v) = ∫((∇ .* (Λf * v)) ⊙ (Λf .* ∇(u)) - (k^2 * ((ε_tag ∘ τ) * v * u)))dΩ + ∫(v * u * w) * dΓ_n # LHS weak form

b_terms = []
for pol in ["X", "Y"]
    for i in 1:N
        src = source_functions[(pol == "X" ? 0 : N) + i]
        RHS(v) = ∫(v * src) * dΩ
        push!(b_terms, RHS)
    end
end

# Solver setup:

print("Solving the linear system(s) of equations...\n")

uh_list = []

for pol in ["X", "Y"]
    for i in 1:N
        println("Solving for dipole $i with polarization $pol...")
        b = b_terms[(pol == "X" ? 0 : N) + i]
        op = AffineFEOperator(a, b, U, V)
        uh = solve(op)
        push!(uh_list, uh)
    end
end

# Plotting fields

print("Outputing results...\n")

if out == true
    for pol in ["X", "Y"]
        for i in 1:N
            idx = (pol == "X" ? 0 : N) + i
            plot_save_field("examples/plots/e_real_field_dip_$(pol)_$(i).png", real(uh_list[idx]), Ω, "seismic")          # Save the real part of electric field
            plot_save_field("examples/plots/e_imag_field_dip_$(pol)_$(i).png", imag(uh_list[idx]), Ω, "seismic")          # Save the imaginary part of electric field
            plot_save_norm("examples/plots/e_norm_field_dip_$(pol)_$(i).png", uh_list[idx], Ω, "inferno", false)           # Save the norm of electric field
        end
    end
    ε_field = CellField(ε_tag∘τ, Ω)
    plot_save_field("examples/plots/perm.png", real(ε_field), Ω, "viridis")             # Save the real part of the permittivity
end

