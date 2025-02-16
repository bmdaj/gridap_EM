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

# We stary by defining the parameters of the problem:

print("Defining model parameters...\n")

λ = 35.0          # Wavelength (arbitrary unit)
k = 2*π/λ        # Wave number 
ε₁ = 4.0        # Relative electric permittivity for the material
ε₀ = 1.0         # Relative electric permittivity for the background
εₛ = [ε₁]                  # List of relative permittivities
tag_list = ["Passive"]   # List of tag names
w_src = -195.0           # Source position
out = true       # Output the results to a file

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
dirichlet_tags = "Edges"

U, V, Ω, dΩ, Γ_n, dΓ_n, _, _ = fea_init(model, order, degree, dirichlet_tags) 

ε_tag, τ = set_tags(model, εₛ, tag_list, ε₀)                                         # We set the permittivity tags for the design and passive regions 

w(x) = -1.0im * k                                                          # absorbing boundary cpndition coefficient
a(u,v) = ∫(  (∇(v))⊙(∇(u))- (k^2*((ε_tag∘τ)*v*u)))dΩ #+ ∫( v*u*w )*dΓ_n # LHS weak form
b(v) = 0.0
a1(u,v) =  ∫(v*u)dΩ 
print("Solving the linear system of equations...\n")

A = get_matrix(AffineFEOperator(a,b,U,V))
#A1 = get_matrix(AffineFEOperator(a1,b,U,V))

nev = 3
λ, ϕ = eigs(A; nev=nev, sigma=-k^2*ε₁)

ϕ_normal = ϕ[:,1]/maximum(abs.(real(ϕ[:,1])))

ϕ_cell = FEFunction(V, ϕ_normal)

fig, ax, plt = GLMakie.plot(Ω, real(ϕ_cell), colormap=:inferno)
Colorbar(fig[1, 2], plt)
save("examples/plots/mode.png", fig)

coords = Gridap.Geometry.get_node_coordinates(Ω)
ϕ_line_values = []
y_values = []
for coord in coords
    x = coord[1]
    if abs(x - w_src) < 1e-10
        push!(ϕ_line_values, ϕ_cell(coord))
        push!(y_values, coord[2])
    end
end
# Reorder y_values in ascending order and order ϕ_line_values accordingly
sorted_indices = sortperm(y_values)
y_values = y_values[sorted_indices]
ϕ_line_values = ϕ_line_values[sorted_indices]



Plots.plot(y_values, real(ϕ_line_values))


#fig, ax, plt = Plots.plot(knots[1], ϕ_mode(knots[1]))
#save("examples/plots/mode_line.png", fig)
