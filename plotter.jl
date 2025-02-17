using Gridap, GridapMakie, GLMakie

function plot_save_field(filename, uh, Ω, cmap)
    fig, ax, plt = GLMakie.plot(Ω, uh, colormap=cmap)
    ax.aspect = DataAspect()
    Colorbar(fig[1, 2], plt)
    save(filename, fig)
end
