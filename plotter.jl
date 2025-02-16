using Gridap, GridapMakie, GLMakie

function plot_e_field(uh, Ω)
    fig, ax, plt = GLMakie.plot(Ω, imag(uh))
    ax.aspect = DataAspect()
    Colorbar(fig[1, 2], plt)
    save("examples/plots/e_field.png", fig)
end

function plot_e_norm(uh, Ω)
    fig, ax, plt = GLMakie.plot(Ω, real(uh * conj(uh)), colormap=:inferno)
    ax.aspect = DataAspect()
    Colorbar(fig[1, 2], plt)
    save("examples/plots/e_norm.png", fig)
end

function plot_perm(ε, Ω)
    fig, ax, plt = GLMakie.plot(Ω, real(ε), colormap=:viridis)
    ax.aspect = DataAspect()
    Colorbar(fig[1, 2], plt)
    save("examples/plots/perm.png", fig)
end