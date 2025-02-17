using Gmsh
import Gmsh: gmsh

function create_waveguide_mesh(w_tot, h_tot, d_tot, w_wg, h_wg, d_wg, w_src, w_dr, h_dr, lc)
    gmsh.initialize()
    gmsh.model.add("waveguide")


    x_off = -w_tot/2
    y_off = -h_tot/2

    # Define points with characteristic length lc and apply offset
    p4 = gmsh.model.geo.addPoint(0 + x_off, 0 + y_off, 0, lc)
    p3 = gmsh.model.geo.addPoint(w_tot + x_off, 0 + y_off, 0, lc)
    p2 = gmsh.model.geo.addPoint(w_tot + x_off, h_tot + y_off, 0, lc)
    p1 = gmsh.model.geo.addPoint(0 + x_off, h_tot + y_off, 0, lc)
    p16 = gmsh.model.geo.addPoint(w_src + x_off, 0 + y_off, 0, lc)
    p13 = gmsh.model.geo.addPoint(w_src + x_off, h_tot + y_off, 0, lc)

    # Define waveguide points and apply offset
    p7 = gmsh.model.geo.addPoint(0 + x_off, (h_tot - h_wg) / 2 + y_off, 0, lc)
    p5 = gmsh.model.geo.addPoint(0 + x_off, (h_tot + h_wg) / 2 + y_off, 0, lc)
    p14 = gmsh.model.geo.addPoint(w_src + x_off, (h_tot + h_wg) / 2 + y_off, 0, lc)
    p15 = gmsh.model.geo.addPoint(w_src + x_off, (h_tot - h_wg) / 2 + y_off, 0, lc)
    p6 = gmsh.model.geo.addPoint(w_wg + x_off, (h_tot + h_wg) / 2 + y_off, 0, lc)
    p8 = gmsh.model.geo.addPoint(w_wg + x_off, (h_tot - h_wg) / 2 + y_off, 0, lc)

    # Define design region points and apply offset
    p12 = gmsh.model.geo.addPoint(w_wg + x_off, (h_tot - h_dr) / 2 + y_off, 0, lc)
    p11 = gmsh.model.geo.addPoint(w_wg + w_dr + x_off, (h_tot - h_dr) / 2 + y_off, 0, lc)
    p10 = gmsh.model.geo.addPoint(w_wg + w_dr + x_off, (h_tot + h_dr) / 2 + y_off, 0, lc)
    p9 = gmsh.model.geo.addPoint(w_wg + x_off, (h_tot + h_dr) / 2 + y_off, 0, lc)


    # Define lines
    l1 = gmsh.model.geo.addLine(p1, p13,1)
    l2 = gmsh.model.geo.addLine(p13, p2,2)
    l3 = gmsh.model.geo.addLine(p2, p3,3)
    l4 = gmsh.model.geo.addLine(p3, p16,4)
    l5 = gmsh.model.geo.addLine(p16, p4,5)
    l6 = gmsh.model.geo.addLine(p4, p7,6)
    l7 = gmsh.model.geo.addLine(p7, p5,7)
    l8 = gmsh.model.geo.addLine(p5, p1,8)
    l9 = gmsh.model.geo.addLine(p5, p14,9)
    l10 = gmsh.model.geo.addLine(p7, p15,10)
    l11 = gmsh.model.geo.addLine(p14, p6,11)
    l12 = gmsh.model.geo.addLine(p15, p8,12)
    l13 = gmsh.model.geo.addLine(p8, p6,13)
    l14 = gmsh.model.geo.addLine(p6, p9,14)
    l15 = gmsh.model.geo.addLine(p9, p10,15)
    l16 = gmsh.model.geo.addLine(p10, p11,16)
    l17 = gmsh.model.geo.addLine(p11, p12,17)
    l18 = gmsh.model.geo.addLine(p12, p8,18)
    l19 = gmsh.model.geo.addLine(p13, p14,19)
    l20 = gmsh.model.geo.addLine(p14, p15,20)
    l21 = gmsh.model.geo.addLine(p15, p16,21)

    # Define curve loops
    cl_tot_1 = gmsh.model.geo.addCurveLoop([l2, l3, l4, -l21, l12, -l18, -l17, -l16, -l15, -l14, -l11, -l19])
    cl_tot_2 = gmsh.model.geo.addCurveLoop([l1, l19, -l9, l8])
    cl_tot_3 = gmsh.model.geo.addCurveLoop([l10, l21, l5, l6])
    cl_wg_1 =  gmsh.model.geo.addCurveLoop([l7, l9, l20, -l10])
    cl_wg_2 =  gmsh.model.geo.addCurveLoop([l11, -l13, -l12, -l20])
    cl_dr = gmsh.model.geo.addCurveLoop([l13, l14, l15, l16, l17, l18])

   
    # Define surfaces for the regions between the waveguide and top/bottom regions
    surf_tot_1 = gmsh.model.geo.addPlaneSurface([cl_tot_1],1)
    surf_tot_2 = gmsh.model.geo.addPlaneSurface([cl_tot_2],2)
    surf_tot_3 = gmsh.model.geo.addPlaneSurface([cl_tot_3],3)
    surf_wg_1 = gmsh.model.geo.addPlaneSurface([cl_wg_1],4)
    surf_wg_2 = gmsh.model.geo.addPlaneSurface([cl_wg_2],5)
    surf_dr = gmsh.model.geo.addPlaneSurface([cl_dr],6)

    surf_list = [surf_tot_1, surf_tot_2, surf_tot_3, surf_wg_1, surf_wg_2, surf_dr]

    # Loop over the surface list and add physical groups
    for surf in surf_list
        
        num_el_1 = 4
        num_el_2 = 8

        fac = d_wg / d_tot
        gmsh.model.geo.extrude([(2, surf)], 0, 0, -d_tot, [num_el_1, num_el_2], [fac, 1.0])
        #ext1 = gmsh.model.geo.extrude([(2, surf)], 0, 0, -d_wg, [num_el])
    
    end

    #gmsh.model.addPhysicalGroup(1, [1,2,3,4,5,6,7,8], 1)
    #gmsh.model.setPhysicalName(1, 1, "Edges")

    #gmsh.model.addPhysicalGroup(1, [19,20,21], 2)
    #gmsh.model.setPhysicalName(1, 2, "Source")

    #gmsh.model.addPhysicalGroup(2, [1,2,3], 3)
    #gmsh.model.setPhysicalName(2, 3, "Air")
    #gmsh.model.addPhysicalGroup(2, [6], 4)
    #gmsh.model.setPhysicalName(2, 4, "Design")
    #gmsh.model.addPhysicalGroup(2, [4,5], 5)
    #gmsh.model.setPhysicalName(2, 5, "Passive")

    # Synchronize and generate mesh
    gmsh.model.geo.synchronize()
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)

    # Save mesh to file
    gmsh.write("waveguide_mesh_3D.msh")
    gmsh.finalize()
end

# Geometry parameters
λ = 35               # Wavelength (arbitrary unit)
w_tot = 400          # Width of the area
h_tot = 200          # Height of the area
d_tot = 100          # Depth of the domain
w_wg = 200           # Width of the waveguide
h_wg = 15            # Height of the waveguide
d_wg = 15            # Depth of the waveguide
w_src  = 80          # Location of the source
w_dr = 75            # Width of the design region
h_dr = 75            # Height of the design region

resol = 5.0         # Number of points per wavelength
lc = λ/resol         # Characteristic length

create_waveguide_mesh(w_tot, h_tot, d_tot, w_wg, h_wg, d_wg, w_src, w_dr, h_dr, lc)