using Gridap, GridapMakie, GLMakie

function plot_e_field(uh, Ω)
    fig, ax, plt = plot(Ω, imag(uh))
    Colorbar(fig[1, 2], plt)
    save("examples/plots/e_field.png", fig)
end

function plot_e_norm(uh, Ω)
    fig, ax, plt = plot(Ω, real(uh * conj(uh)), colormap=:inferno)
    Colorbar(fig[1, 2], plt)
    save("examples/plots/e_norm.png", fig)
end

function plot_perm(ε, Ω)
    fig, ax, plt = plot(Ω, real(ε), colormap=:viridis)
    Colorbar(fig[1, 2], plt)
    save("examples/plots/perm.png", fig)
end