using GridapGmsh
using Gridap.Fields

include("../msh/mesh_metalens.jl")
include("../fea.jl")
include("../topopt.jl")
include("../pml.jl")
include("../plotter.jl")

# Geometry parameters
λ = 35           # Wavelength (arbitrary unit)
L = 200          # Width of the area
H = 200          # Height of the area
hs = -2.0 * 35  # y-position of the source (plane wave)
h = 2.5 * 35    # displacement of metalens with respect to end
hsub = 5        # subtrate thickness
hd = 15          # thickness of design region
d_pml = 35       # Thickness of the PML

resol = 5.0      # Number of points per wavelength
lc = λ/resol      # Characteristic length

#MeshGenerator(L,H, hs, h, hsub, hd, d_pml,lc)

# We stary by defining the parameters of the problem:

print("Defining model parameters...\n")

k = 2*π/λ                          # Wave number 
ε₁ = 4.0                           # Relative electric permittivity for the material
ε₀ = 1.0                           # Relative electric permittivity for the background
εₛ = [ε₁]                          # List of relative permittivities
tag_list = ["Passive"]             # List of tag names
out = true                         # Output the results to a file

d_pml = λ                          # Thickness of the PML
w_tot = L  + 2 * d_pml
h_tot = H  + 2 * d_pml

# We define the optimization parameters:

print("Defining optimization parameters...\n")

r_f = (λ/6) /sqrt(3)               # filter radius for PDE filter
β = 5.0                            # threshold sharpness for smoothed Heaviside threshold
η = 0.5                            # threshold value for smoothed Heaviside threshold
ξ₀ = 0.5                           # Homogeneous initial guess for the design variables

# We import the geometry of the problem from the file `geometry.msh`:

print("Importing mesh...\n")

model = GmshDiscreteModel("geometry.msh")

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
source_tags =  "Source"
design_tags = "Design"

U, V, Ω, dΩ, Ω_d, dΩ_d, Γ_n, dΓ_n, Γ_s, dΓ_s = fea_init_topopt(model, order, degree, dirichlet_tags, neumann_tags, design_tags, source_tags) 

# We set up the design parameters and filter and threshold the design variables:

design_params = design_variables_init()
ξ_0_vec =  ξ₀ * ones(design_params.np)
ξ_f_vec = Filter(ξ_0_vec; r_f, dΩ_d)
ξ_fh = FEFunction(design_params.Pf, ξ_f_vec)
ξ_th = (ξ_f -> Threshold(ξ_f; β, η)) ∘ ξ_fh

# We set the permittivity tags for non-design region: 

ε_tag, τ = set_tags(model, εₛ, tag_list, ε₀)

# We set the PML function:

Λf = pml(w_tot, h_tot, k, d_pml)         

# We set the weak form of the PDE:

w(x) = -1.0im * k                                                                                          # absorbing boundary condition coefficient
a₀(u, v) = ∫((∇ .* (Λf * v)) ⊙ (Λf .* ∇(u)) - (k^2 * ((ε_tag ∘ τ) * v * u)))dΩ + ∫(v * u * w) * dΓ_n       # LHS weak form for base problem
a_ξ(u,v,ξ_th) = ∫(∇(v)⊙ (∇(u)) - (k^2 * (((ξ -> ξ_inter(ξ, ε₀, ε₁)) ∘ ξ_th) * v * u)))dΩ_d      # LHS weak form in design domain

# Matrix assembly:

print("Assembling the matrices and vectors...\n")

function MatrixA(ξ_th; U, V)

    A_mat = assemble_matrix(U, V) do u, v
        a₀(u, v) + a_ξ(u, v, ξ_th)
    end
    
    return lu(A_mat)
    
end

start_time = time_ns()  # Get start time in nanoseconds

A_mat = MatrixA(ξ_th; U, V)
b_vec = assemble_vector(v->(∫(v)dΓ_s), V)

# Solver setup:

print("Solving the linear system of equations...\n")

Ez_vec = A_mat \ b_vec
Ez = FEFunction(U, Ez_vec)

end_time = time_ns()  # Get end time

elapsed_time = (end_time - start_time) / 1e9  # Convert to seconds
println("Elapsed time: $elapsed_time seconds")


using Plots

times = [5.78, 6.06, 6.74, 9.39, 23.08, 113.25]
size = [1/5, 1/10, 1/25, 1/50, 1/100, 1/200]
DOFs = [2025, 7392, 88123, 350776, 1387392, 2759031]

# Plotting the results
p1 = Plots.scatter(size, times, xaxis=(:log, :flip), xlabel="Size (λ) ", ylabel="Time (s)", label="Time", legend=:topleft)
p2 = Plots.scatter(DOFs, times, xaxis=:log, xlabel="DOFs", ylabel="Time (s)", label="Time", legend=:topleft)

Plots.plot(p1, p2, layout=(2, 1))

savefig("examples/plots/scaling_single_core_metalens.png")


