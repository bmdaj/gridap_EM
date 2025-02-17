using GridapGmsh
using Gridap.Fields
using Interpolations
using Arpack

include("../fea.jl")
include("../plotter.jl")
include("../pml.jl")
# We stary by defining the parameters of the problem:

print("Defining model parameters...\n")

λ = 35.0                           # Wavelength (arbitrary unit)
k = 2 * π / λ                          # Wave number 
ε₁ = 4.0                           # Relative electric permittivity for the material
ε₀ = 1.0                           # Relative electric permittivity for the background
εₛ = [0.5 * ε₁, ε₁]                  # List of relative permittivities
tag_list = ["Design", "Passive"]   # List of tag names
out = true                         # Output the results to a file

w_tot = 400.0
h_tot = 200.0
d_pml = λ                          # Thickness of the PML


# We import the input mode to set it as a source:

print("Importing input mode...\n")

data = readdlm("examples/data/mode_line.txt", Float64)
y_values = data[:, 1]
ϕ_line_values_real = data[:, 2]
ϕ_line_values_imag = data[:, 3]
ϕ_line_values = ϕ_line_values_real + 1im*ϕ_line_values_imag

ϕ_interp = Interpolations.interpolate((y_values,), ϕ_line_values, Gridded(Linear()))
ϕ_mode(x) = ϕ_interp(x[2])

# We import the geometry of the problem from the file `geometry.msh`:

print("Importing mesh...\n")

model = GmshDiscreteModel("waveguide_mesh.msh")

print("Outputing mesh...\n")

if out == true
    writevtk(model, "geom/geometry")
end

# FEA setup:

print("Setting up FEA...\n")

order = 1
degree = 2
dirichlet_tags = "None"
neumann_tags = "Edges"
source_tags = "Source"

U, V, Ω, dΩ, Γ_n, dΓ_n, Γ_s, dΓ_s = fea_init(model, order, degree, dirichlet_tags, neumann_tags, source_tags)

ε_tag, τ = set_tags(model, εₛ, tag_list, ε₀)                                         # We set the permittivity tags for the design and passive regions 

Λf = pml(w_tot, h_tot, k, d_pml)                                                               # We set the PML function

w(x) = -1.0im * k                                                          # absorbing boundary cpndition coefficient
a(u, v) = ∫((∇ .* (Λf * v)) ⊙ (Λf .* ∇(u)) - (k^2 * ((ε_tag ∘ τ) * v * u)))dΩ + ∫(v * u * w) * dΓ_n # LHS weak form

b(v) = ∫(ϕ_mode * v) * dΓ_s                                                           # RHS weak form

# Solver setup:

print("Solving the linear system of equations...\n")

op = AffineFEOperator(a, b, U, V)
uh = solve(op)

# Plotting fields

print("Outputing results...\n")

if out == true

    plot_save_field("examples/plots/e_real_field.png", real(uh), Ω, "seismic")          # Save the real part of electric field
    plot_save_field("examples/plots/e_imag_field.png", imag(uh), Ω, "seismic")          # Save the imaginary part of electric field
    plot_save_field("examples/plots/e_norm_field.png", real(uh*conj(uh)), Ω, "inferno") # Save the imaginary part of electric field

    ε_field = CellField(ε_tag∘τ, Ω)

    plot_save_field("examples/plots/perm.png", real(ε_field), Ω, "viridis")             # Save the real part of the permittivity

end

