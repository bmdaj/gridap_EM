using Gmsh
import Gmsh: gmsh
function MeshGenerator(L,H, hs, h, hsub, hd, d_pml,lc)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.clear()
    gmsh.model.add("geometry")

    # Add points
    gmsh.model.geo.addPoint(-L/2-d_pml, -H/2-d_pml, 0, lc, 1)
    gmsh.model.geo.addPoint( L/2+d_pml, -H/2-d_pml, 0, lc, 2)
    gmsh.model.geo.addPoint( L/2+d_pml,  hs       , 0, lc, 3)
    gmsh.model.geo.addPoint(-L/2-d_pml,  hs       , 0, lc, 4)
    gmsh.model.geo.addPoint( L/2+d_pml,  H/2+d_pml, 0, lc, 5)
    gmsh.model.geo.addPoint(-L/2-d_pml,  H/2+d_pml, 0, lc, 6)
    gmsh.model.geo.addPoint(-L/2      ,  -H/2 + hsub + h, 0, lc, 7)
    gmsh.model.geo.addPoint(-L/2      ,  -H/2 + hsub + hd + h, 0, lc, 8)
    gmsh.model.geo.addPoint( L/2      ,  -H/2 + hsub + h , 0, lc, 9)
    gmsh.model.geo.addPoint( L/2      ,  -H/2 + hsub + hd + h, 0, lc, 10)
    gmsh.model.geo.addPoint( -L/2     ,  -H/2 + h, 0, lc, 11)
    gmsh.model.geo.addPoint( L/2      ,  -H/2 + h, 0, lc, 12)

    # Add lines
    gmsh.model.geo.addLine( 1,  2,  1)
    gmsh.model.geo.addLine( 2,  3,  2)
    gmsh.model.geo.addLine( 3,  4,  3)
    gmsh.model.geo.addLine( 1,  4,  4)
    gmsh.model.geo.addLine( 3,  5,  5)
    gmsh.model.geo.addLine( 5,  6,  6)
    gmsh.model.geo.addLine( 4,  6,  7)

    gmsh.model.geo.addLine( 7,  8,  8)
    gmsh.model.geo.addLine( 11,  7,  9)
    gmsh.model.geo.addLine( 11,  12,  10)
    gmsh.model.geo.addLine( 8,  10,  11)
    gmsh.model.geo.addLine( 7,  9,  12)
    gmsh.model.geo.addLine( 12,  9,  13)
    gmsh.model.geo.addLine( 9,  10,  14)

    # Construct curve loops and surfaces 
    gmsh.model.geo.addCurveLoop([1, 2, 3, -4], 1)
    gmsh.model.geo.addCurveLoop([5, 6,-7, -3], 2)
    gmsh.model.geo.addCurveLoop([-8, -11, 14, 12], 3)
    gmsh.model.geo.addCurveLoop([-9, -12, 13, 10], 4)

    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addPlaneSurface([2,3,4], 2)
    gmsh.model.geo.addPlaneSurface([3], 3)
    gmsh.model.geo.addPlaneSurface([4], 4)
    # Physical groups
    gmsh.model.addPhysicalGroup(1, [1,2,4,5,6,7], 1)
    #gmsh.model.addPhysicalGroup(1, [2,4,5,6,7], 1)
    gmsh.model.setPhysicalName(1, 1, "Edges")
    gmsh.model.addPhysicalGroup(1, [3], 2)
    #gmsh.model.addPhysicalGroup(1, [1], 2)
    gmsh.model.setPhysicalName(1, 2, "Source")

    gmsh.model.addPhysicalGroup(2, [1,2], 3)
    gmsh.model.setPhysicalName(2, 3, "Air")
    gmsh.model.addPhysicalGroup(2, [3], 4)
    gmsh.model.setPhysicalName(2, 4, "Design")
    gmsh.model.addPhysicalGroup(2, [4], 5)
    gmsh.model.setPhysicalName(2, 5, "Passive")

    gmsh.model.geo.synchronize()

    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)
    # ... and save it to disk
    gmsh.write("geometry.msh")
    gmsh.finalize()
end

# Geometry parameters
λ = 35           # Wavelength (arbitrary unit)
L = 200          # Width of the area
H = 200          # Height of the area
hs = -2.0 * 35  # y-position of the source (plane wave)
h = 2.5 * 35    # displacement of metalens with respect to end
hsub = 5        # subtrate thickness
hd = 15          # thickness of design region
d_pml = 35       # Thickness of the PML

resol = 25.0      # Number of points per wavelength
lc = λ/resol      # Characteristic length

MeshGenerator(L,H, hs, h, hsub, hd, d_pml,lc)