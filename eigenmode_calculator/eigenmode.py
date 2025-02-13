import meep as mp
from meep import mpb

def eigenmode_1D(eps_back, eps_wg, thk_wg, thk_dom, res, k_point, bands, pol="TE"):


    geometry_lattice = mp.Lattice(size=(0, thk_dom))

    default_material = mp.Medium(epsilon=eps_back)

    geometry = [mp.Block(center=(0,0), # center of computational cell
                     size=(mp.inf, thk_wg, mp.inf),
                     material=mp.Medium(epsilon=eps_wg))]


    k_points =[mp.Vector3(k_point)]

    ms = mpb.ModeSolver(num_bands=bands,
                    k_points=k_points,
                    geometry=geometry,
                    geometry_lattice=geometry_lattice,
                    resolution=res,
                    default_material=default_material)
    
    if pol == "TE":
        ms.run_te(mpb.fix_efield_phase)
    elif pol == "TM":
        raise ValueError("TM polarization not implemented yet")
    
    return ms
