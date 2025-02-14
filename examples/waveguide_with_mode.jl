using GridapGmsh
using Gridap.Fields
using Interpolations
using Arpack

include("../fea.jl")
include("../plotter.jl")

# We stary by defining the parameters of the problem:

print("Defining model parameters for the mode calculation...\n")

λ = 35.0          # Wavelength (arbitrary unit)
k = 2*π/λ        # Wave number 
ε₁ = 4.0        # Relative electric permittivity for the material
ε₀ = 1.0         # Relative electric permittivity for the background
εₛ = [ε₁]                  # List of relative permittivities
tag_list = ["Passive"]   # List of tag names
out = true       # Output the results to a file

# We import the geometry of the problem from the file `geometry.msh`:

print("Importing mode calculation mesh...\n")

model = GmshDiscreteModel("waveguide_mode_mesh.msh")

print("Outputing mode calculation mesh...\n")

if out == true
    writevtk(model,"geom/geometry")
end

# FEA setup:

print("Setting up mode FEA...\n")

order = 1 
degree = 2
dirichlet_tags = "Edges"

U, V, Ω, dΩ, Γ_n, dΓ_n, _, _ = fea_init(model, order, degree, dirichlet_tags) 

ε_tag, τ = set_tags(model, εₛ, tag_list, ε₀)                                         # We set the permittivity tags for the design and passive regions 

w(x) = -1.0im * k                                                          # absorbing boundary cpndition coefficient
a(u,v) = ∫(  (∇(v))⊙(∇(u))- (k^2*((ε_tag∘τ)*v*u)))dΩ 
b(v) = 0.0
a1(u,v) =  ∫(v*u)dΩ 
print("Solving the eigenvalue problem...\n")

A = get_matrix(AffineFEOperator(a,b,U,V))

nev = 3
λ, ϕ = eigs(A; nev=nev, sigma=-k^2*ε₁)

ϕ_cell = FEFunction(V, ϕ[:,1])
coords = Gridap.Geometry.get_node_coordinates(Ω)
ϕ_line_values = []
y_values = []
for coord in coords
    x = coord[1]
    if abs(x - 5.0) < 1e-10
        push!(ϕ_line_values, ϕ_cell(coord))
        push!(y_values, coord[2])
    end
end
# Reorder y_values in ascending order and order ϕ_line_values accordingly
sorted_indices = sortperm(y_values)
y_values = y_values[sorted_indices]
ϕ_line_values = ϕ_line_values[sorted_indices]
ϕ_interp = Interpolations.interpolate((y_values,), ϕ_line_values, Gridded(Linear()))
ϕ_mode(x) = ϕ_interp(x[2])


#ϕ_mode(x)= ϕ[x[1], x[2]]

# We stary by defining the parameters of the problem:

print("Defining model parameters...\n")

λ = 35.0                           # Wavelength (arbitrary unit)
k = 2*π/λ                          # Wave number 
ε₁ = 4.0                           # Relative electric permittivity for the material
ε₀ = 1.0                           # Relative electric permittivity for the background
εₛ = [0.5*ε₁, ε₁]                  # List of relative permittivities
tag_list = ["Design", "Passive"]   # List of tag names
out = true                         # Output the results to a file

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
source_tags = "Source"

U, V, Ω, dΩ, Γ_n, dΓ_n, Γ_s, dΓ_s = fea_init(model, order, degree, dirichlet_tags, neumann_tags, source_tags) 

ε_tag, τ = set_tags(model, εₛ, tag_list, ε₀)                                         # We set the permittivity tags for the design and passive regions 

w(x) = -1.0im * k                                                          # absorbing boundary cpndition coefficient
a(u,v) = ∫(  (∇(v))⊙(∇(u)) - (k^2*((ε_tag∘τ)*v*u))  )dΩ  + ∫( v*u*w )*dΓ_n # LHS weak form

b(v) = ∫(ϕ_mode*v)*dΓ_s                                                           # RHS weak form

# Solver setup:

print("Solving the linear system of equations...\n")

op = AffineFEOperator(a,b,U,V)
uh = solve(op)

# Plotting fields

print("Outputing results...\n")

if out == true

    plot_e_field(uh,Ω)
    plot_e_norm(uh,Ω)

    ε_field = CellField(ε_tag∘τ, Ω)
    plot_perm(ε_field,Ω) 
    
end

