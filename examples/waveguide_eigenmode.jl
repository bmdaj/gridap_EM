using GridapGmsh
using Gridap.Fields
using Gridap
using Arpack 
using DelimitedFiles
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.Fields
using Plots

include("../fea.jl")
include("../plotter.jl")

# We start by defining the parameters of the problem:

print("Defining model parameters...\n")

λ = 35.0                    # Wavelength (arbitrary unit)
k = 2*π/λ                   # Wave number 
ε₁ = 4.0                    # Relative electric permittivity for the material
ε₀ = 1.0                    # Relative electric permittivity for the background
εₛ = [ε₁]                   # List of relative permittivities
tag_list = ["Passive"]      # List of tag names
out = true                  # Output the results to a file

# We the define the parameters for the mode solver:

print("Defining mode solver parameters...\n")

w_src = -195.0              # Source position
nev = 3
n_mode = 1

# We import the geometry of the problem from the file `geometry.msh`:

print("Importing mesh...\n")

model = GmshDiscreteModel("waveguide_mode_mesh.msh")

print("Outputing mesh...\n")

if out == true
    writevtk(model,"geom/geometry")
end

# FEA setup:

print("Setting up FEA...\n")

order = 1 
degree = 2
dirichlet_tags = "Edges"                                 # We set PEC (dirichlet) boundary conditions on the ends of the simulation domain

U, V, Ω, dΩ, Γ_n, dΓ_n, _, _ = fea_init(model, order, degree, dirichlet_tags) 

ε_tag, τ = set_tags(model, εₛ, tag_list, ε₀)             # We set the permittivity tags for the design and passive regions 

a(u,v) = ∫(  (∇(v))⊙(∇(u))- (k^2*((ε_tag∘τ)*v*u)))dΩ     # Note that there is no absorbing boundary conditions
b(v) = 0.0                                               # There is no source for the eigenvalue problem
print("Solving the linear system of equations...\n")

A = get_matrix(AffineFEOperator(a,b,U,V))
λ, ϕ = eigs(A; nev=nev, sigma=-k^2*ε₁)                   # We solve around the propagation constant of the material: sigma=-k^2*ε₁

# Postprocessing and outputing the mode:

print("Postprocessing and outputing results...\n")

ϕ_normal = ϕ[:,n_mode]/maximum(abs.(real(ϕ[:,n_mode])))  # We normalize the mode
ϕ_cell = FEFunction(V, ϕ_normal)

if out 
    plot_save_field("examples/plots/e_real_wg_mode.png", real(ϕ_cell), Ω, "seismic")              # Save the real part of electric field
    plot_save_field("examples/plots/e_imag_wg_mode.png", imag(ϕ_cell), Ω, "seismic")              # Save the imaginary part of electric field
    plot_save_field("examples/plots/e_norm_wg_mode.png", real(ϕ_cell*conj(ϕ_cell)), Ω, "inferno") # Save the imaginary part of electric field
end

# We find the values of the field lying on the source line 

coords = Gridap.Geometry.get_node_coordinates(Ω)
ϕ_line_values = ComplexF64[]
y_values = ComplexF64[]

for coord in coords
    if isapprox(coord[1], w_src; atol=1e-10)
        push!(ϕ_line_values, ϕ_cell(coord))
        push!(y_values, coord[2])
    end
end

# Reorder y_values in ascending order and order ϕ_line_values accordingly

sorted_indices = sortperm(real(y_values))
y_values = y_values[sorted_indices]
ϕ_line_values = ϕ_line_values[sorted_indices]

if out
    writedlm("examples/data/mode_line.txt", hcat(real(y_values), real(ϕ_line_values), imag(ϕ_line_values)))
end

Plots.plot(real(y_values), real(ϕ_line_values))