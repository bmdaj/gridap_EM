using GridapGmsh
using Gridap.Fields
using Gridap
using Arpack 
using DelimitedFiles
using Gridap.FESpaces

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

ϕ_cell = FEFunction(V, ϕ[:,1])


fig, ax, plt = plot(Ω, real(ϕ_cell), colormap=:viridis)
save("examples/plots/mode.png", fig)

xe = get_cell_coordinates(Ω)