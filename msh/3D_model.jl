using Gmsh
import Gmsh: gmsh

function MeshGenerator(L, lc)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
    #gmsh.option.setNumber("Mesh.RecombineAll", 1)  # Recombine all elements

    gmsh.clear()
    gmsh.model.add("cube")

    # Add points
    gmsh.model.geo.addPoint(-L/2, -L/2, -L/2, lc, 1)
    gmsh.model.geo.addPoint(L/2, -L/2, -L/2, lc, 2)
    gmsh.model.geo.addPoint(L/2, L/2, -L/2, lc, 3)
    gmsh.model.geo.addPoint(-L/2, L/2, -L/2, lc, 4)
    gmsh.model.geo.addPoint(-L/2, -L/2, L/2, lc, 5)
    gmsh.model.geo.addPoint(L/2, -L/2, L/2, lc, 6)
    gmsh.model.geo.addPoint(L/2, L/2, L/2, lc, 7)
    gmsh.model.geo.addPoint(-L/2, L/2, L/2, lc, 8)
    gmsh.model.geo.addPoint(0, 0, 0, lc, 9)  # Center point

    # Add lines
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addLine(5, 6, 5)
    gmsh.model.geo.addLine(6, 7, 6)
    gmsh.model.geo.addLine(7, 8, 7)
    gmsh.model.geo.addLine(8, 5, 8)
    gmsh.model.geo.addLine(1, 5, 9)
    gmsh.model.geo.addLine(2, 6, 10)
    gmsh.model.geo.addLine(3, 7, 11)
    gmsh.model.geo.addLine(4, 8, 12)

    # Construct curve loops and surfaces
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 2)
    gmsh.model.geo.addCurveLoop([1, 10, -5, -9], 3)
    gmsh.model.geo.addCurveLoop([2, 11, -6, -10], 4)
    gmsh.model.geo.addCurveLoop([3, 12, -7, -11], 5)
    gmsh.model.geo.addCurveLoop([4, 9, -8, -12], 6)

    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.addPlaneSurface([3], 3)
    gmsh.model.geo.addPlaneSurface([4], 4)
    gmsh.model.geo.addPlaneSurface([5], 5)
    gmsh.model.geo.addPlaneSurface([6], 6)

    # Create volume
    gmsh.model.geo.addSurfaceLoop([1, 2, 3, 4, 5, 6], 1)
    gmsh.model.geo.addVolume([1], 1)

    # Physical groups
    gmsh.model.addPhysicalGroup(0, [9], 1)
    gmsh.model.setPhysicalName(0, 1, "Center Point")
    gmsh.model.addPhysicalGroup(3, [1], 2)
    gmsh.model.setPhysicalName(3, 2, "Cube Volume")

    gmsh.model.geo.synchronize()

    # Generate a 3D mesh
    gmsh.model.mesh.generate(3)
    # Save it to disk
    gmsh.write("cube.msh")
    gmsh.finalize()
end

# Geometry parameters
L = 200          # Side length of the cube
resol = 20.0     # Number of points per unit length
lc = L/resol     # Characteristic length

MeshGenerator(L, lc)
