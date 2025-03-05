using Gridap, GridapMakie, GLMakie, Gridap.Fields

function plot_save_field(filename, uh, Ω, cmap)
    fig, ax, plt = GLMakie.plot(Ω, uh, colormap=cmap)
    ax.aspect = DataAspect()
    Colorbar(fig[1, 2], plt)
    save(filename, fig)
end

function plot_save_norm(filename, uh, Ω, cmap, logscale=false)
    if logscale
        uh_cell = get_free_dof_values(uh)
        uh_cell = log10.(uh_cell.*conj(uh_cell))
        uh = FEFunction(V, uh_cell)
    else
        uh = uh * conj(uh)
    end
    fig, ax, plt = GLMakie.plot(Ω, real(uh), colormap=cmap)
    ax.aspect = DataAspect()
    Colorbar(fig[1, 2], plt)
    save(filename, fig)
end