using Gmsh
import Gmsh: gmsh

function create_waveguide_mesh(w_tot, h_tot, w_wg, h_wg, w_src, w_dr, h_dr, lc)
    gmsh.initialize()
    gmsh.model.add("waveguide")

    x_off = -w_tot/2
    y_off = -h_tot/2

    # Define points with characteristic length lc and apply offset
    p4 = gmsh.model.geo.addPoint(0 + x_off, 0 + y_off, 0, lc)
    p1 = gmsh.model.geo.addPoint(0 + x_off, h_tot + y_off, 0, lc)
    p16 = gmsh.model.geo.addPoint(w_src + x_off, 0 + y_off, 0, lc)
    p13 = gmsh.model.geo.addPoint(w_src + x_off, h_tot + y_off, 0, lc)

    # Define waveguide points and apply offset
    p7 = gmsh.model.geo.addPoint(0 + x_off, (h_tot - h_wg) / 2 + y_off, 0, lc)
    p5 = gmsh.model.geo.addPoint(0 + x_off, (h_tot + h_wg) / 2 + y_off, 0, lc)
    p14 = gmsh.model.geo.addPoint(w_src + x_off, (h_tot + h_wg) / 2 + y_off, 0, lc)
    p15 = gmsh.model.geo.addPoint(w_src + x_off, (h_tot - h_wg) / 2 + y_off, 0, lc)

    # Define lines
    l1 = gmsh.model.geo.addLine(p1, p13,1)
    l5 = gmsh.model.geo.addLine(p16, p4,5)
    l6 = gmsh.model.geo.addLine(p4, p7,6)
    l7 = gmsh.model.geo.addLine(p7, p5,7)
    l8 = gmsh.model.geo.addLine(p5, p1,8)
    l9 = gmsh.model.geo.addLine(p5, p14,9)
    l10 = gmsh.model.geo.addLine(p7, p15,10)
    l19 = gmsh.model.geo.addLine(p13, p14,19)
    l20 = gmsh.model.geo.addLine(p14, p15,20)
    l21 = gmsh.model.geo.addLine(p15, p16,21)

    # Define curve loops
    cl_tot_2 = gmsh.model.geo.addCurveLoop([l1, l19, -l9, l8])
    cl_tot_3 = gmsh.model.geo.addCurveLoop([l10, l21, l5, l6])
    cl_wg_1 =  gmsh.model.geo.addCurveLoop([l7, l9, l20, -l10])

   
    # Define surfaces for the regions between the waveguide and top/bottom regions
    #surf_tot_1 = gmsh.model.geo.addPlaneSurface([cl_tot_1],1)
    surf_tot_2 = gmsh.model.geo.addPlaneSurface([cl_tot_2],2)
    surf_tot_3 = gmsh.model.geo.addPlaneSurface([cl_tot_3],3)
    surf_wg_1 = gmsh.model.geo.addPlaneSurface([cl_wg_1],4)

    gmsh.model.addPhysicalGroup(1, [1,5], 1)
    gmsh.model.setPhysicalName(1, 1, "Edges")

    gmsh.model.addPhysicalGroup(2, [2,3], 3)
    gmsh.model.setPhysicalName(2, 3, "Air")
    gmsh.model.addPhysicalGroup(2, [4], 5)
    gmsh.model.setPhysicalName(2, 5, "Passive")

    # Synchronize and generate mesh
    gmsh.model.geo.synchronize()
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)

    # Save mesh to file
    gmsh.write("waveguide_mode_mesh.msh")
    gmsh.finalize()
end

# Geometry parameters
λ = 35               # Wavelength (arbitrary unit)
w_tot = 400          # Width of the area
h_tot = 200          # Height of the area
w_wg = 200            # Width of the waveguide
h_wg = 15            # Height of the waveguide
w_src  = 5          # Location of the source
w_dr = 75               # Width of the design region
h_dr = 75               # Height of the design region

resol = 100.0      # Number of points per wavelength
lc = λ/resol      # Characteristic length

create_waveguide_mesh(w_tot, h_tot, w_wg, h_wg, w_src, w_dr, h_dr, lc)