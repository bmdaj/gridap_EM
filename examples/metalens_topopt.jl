using GridapGmsh
using Gridap.Fields
using NLopt
using ChainRulesCore, Zygote
import ChainRulesCore: rrule

NO_FIELDS = ZeroTangent()

include("../fea.jl")
include("../topopt.jl")
include("../sensitivity.jl")
include("../plotter.jl")

# We stary by defining the parameters of the problem:

print("Defining model parameters...\n")

λ = 35.0                           # Wavelength (arbitrary unit)
k = 2*π/λ                          # Wave number 
ε₁ = 4.0                           # Relative electric permittivity for the material
ε₀ = 1.0                           # Relative electric permittivity for the background
εₛ = [ε₁]                          # List of relative permittivities
tag_list = ["Passive"]             # List of tag names
out = true                         # Output the results to a file

d_pml = λ                          # Thickness of the PML
w_tot = 200.0 + 2 * d_pml
h_tot = 200.0 + 2 * d_pml

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

A_mat = MatrixA(ξ_th; U, V)
b_vec = assemble_vector(v->(∫(v)dΓ_s), V)

# Solver setup:

print("Solving the linear system of equations...\n")

Ez_vec = A_mat \ b_vec
Ez = FEFunction(U, u_vec)

# Plotting fields

print("Outputing results...\n")

if out == true

    plot_save_field("examples/plots/init_e_real_field.png", real(Ez), Ω, "seismic")          # Save the real part of electric field
    plot_save_field("examples/plots/init_e_imag_field.png", imag(Ez), Ω, "seismic")          # Save the imaginary part of electric field
    plot_save_field("examples/plots/init_e_norm_field.png", real(Ez*conj(Ez)), Ω, "inferno") # Save the imaginary part of electric field

    ε_field = CellField(ε_tag∘τ, Ω) + (ξ -> ξ_inter(ξ, ε₀, ε₁))∘ξ_th

    plot_save_field("examples/plots/init_perm.png", real(ε_field), Ω, "viridis")             # Save the real part of the permittivity
    
end

# Objective: optimize the field intensity at the desired position

println("Setting up the FOM and the sensitivities...")

function MatrixFOM(U, V, dΩ)

    x0 = VectorValue(0,50)  # Position of the field to be optimized
    δ = 1

    return assemble_matrix(U, V) do u, v
        ∫((x->(1/(2*π)*exp(-norm(x - x0)^2 / 2 / δ^2))) * (∇(u) ⋅ ∇(v)) )dΩ
    end
end

function FOM_eval(ξ_f_vec; U, V, β, η, dΩ, dΓ_s, design_params)

    ξ_fh = FEFunction(design_params.Pf, ξ_f_vec)
    ξ_th = (ξ_f -> Threshold(ξ_f; β, η)) ∘ ξ_fh
    A_mat = MatrixA(ξ_th; U, V)
    b_vec = assemble_vector(v->(∫(v)dΓ_s), V)
    Ez_vec = A_mat \ b_vec

    O_mat = MatrixFOM(U, V, dΩ)
    real(Ez_vec' * O_mat * Ez_vec)
end

function rrule(::typeof(FOM_eval), ξ_f_vec; U, V, β, η, dΩ, dΓ_s, design_params)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * dFOM_dξf(ξ_f_vec; U, V, β, η, design_params)
    end
    FOM_eval(ξ_f_vec; U, V, β, η, dΩ, dΓ_s, design_params), U_pullback
end

function dFOM_dξf(ξ_f_vec; U, V, β, η, design_params)
    
    ξ_fh = FEFunction(design_params.Pf, ξ_f_vec)
    ξ_th = (ξ_f -> Threshold(ξ_f; β, η)) ∘ ξ_fh
    A_mat = MatrixA(ξ_th; U, V)
    b_vec = assemble_vector(v->(∫(v)dΓ_s), V)
    Ez_vec = A_mat \ b_vec

    O_mat = MatrixFOM(U, V, dΩ)

    Ez = FEFunction(U, Ez_vec)
    λ_adj_vec =  A_mat' \ (O_mat * Ez_vec)
    λ_adj_conjh = FEFunction(U, conj(λ_adj_vec))

    l_temp(dξ) = ∫(real(-2 * dA_dξf(Ez, λ_adj_conjh, ξ_fh; β, η)) * dξ)dΩ_d
    dfom_dξf = assemble_vector(l_temp, design_params.Pf)

    return dfom_dξf
end

ξ_test = rand(design_params.np)
δξ₀ = 1e-8
δξ = δξ₀*rand(design_params.np)
grad = zeros(design_params.np)

FOM₀ = dFOM_dξ(ξ_test, grad; U, V, r_f, β, η, dΩ, dΩ_d, dΓ_s, design_params)
FOM₁ = dFOM_dξ(ξ_test+δξ, []; U, V, r_f, β, η, dΩ, dΩ_d, dΓ_s, design_params)

println("Finite difference check...")
println(FOM₁-FOM₀)
println(grad'*δξ)

println("Setting up the optimization...")

function dFOM_dξ_optimize(ξ₀; r_f, β, η, TOL = 1e-4, MAX_ITER = 500) # NEED TO ADD THE CORRECT ARGUMENTS HERE
    ##################### Optimize #################
    opt = Opt(:LD_MMA, design_params.np)
    opt.lower_bounds = 0
    opt.upper_bounds = 1
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (ξ, grad) -> dFOM_dξ(ξ, grad; U, V, r_f, β, η, dΩ, dΩ_d, dΓ_s, design_params)

    (FOM_opt, ξ_opt, ret) = optimize(opt, ξ₀)
    @show numevals = opt.numevals # the number of function evaluations
    return FOM_opt, ξ_opt
end

# STILL NEED TO ADD THE CONTINUATIOn SCHEME
FOM_opt = 0
ξ_opt = ξ_0_vec
TOL = 1e-8
MAX_ITER = 100

FOM_opt, ξ_opt = dFOM_dξ_optimize(ξ_opt; r_f, β, η, TOL, MAX_ITER)

println("Obtaining optimized solution...")

ξ_f_vec = ξf_ξ(ξ_opt; r_f, dΩ_d)
ξ_fh = FEFunction(design_params.Pf, ξ_f_vec)
ξ_th = (ξ_f -> Threshold(ξ_f; β, η)) ∘ ξ_fh
A_mat = MatrixA(ξ_th; U, V)
b_vec = assemble_vector(v->(∫(v)dΓ_s), V)
Ez_vec = A_mat \ b_vec
Ez = FEFunction(U, Ez_vec)

if out == true

    plot_save_field("examples/plots/opt_e_real_field.png", real(Ez), Ω, "seismic")          # Save the real part of electric field
    plot_save_field("examples/plots/opt_e_imag_field.png", imag(Ez), Ω, "seismic")          # Save the imaginary part of electric field
    plot_save_field("examples/plots/opt_e_norm_field.png", real(Ez*conj(Ez)), Ω, "inferno") # Save the imaginary part of electric field

    ε_field = CellField(ε_tag∘τ, Ω) + (ξ -> ξ_inter(ξ, ε₀, ε₁))∘ξ_th

    plot_save_field("examples/plots/opt_perm.png", real(ε_field), Ω, "viridis")             # Save the real part of the permittivity
    
end