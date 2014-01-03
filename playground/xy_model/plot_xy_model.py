import numpy as np
import matplotlib.pyplot as plt


def main():
    from xy_model_system import XYModlelSystem
    fname = "gmin.coords"
    fname = "goldstone_mode.coords"
    fname = "vortex2.coords"
    reduced_coords = np.genfromtxt(fname)
    print (len(reduced_coords))
    print np.sqrt(len(reduced_coords))
    dim = [24, 24]
    print "dim", dim
    
    system = XYModlelSystem(dim, phi_disorder=0)
    coords = system.pot.coords_converter.get_full_coords(reduced_coords)
    
    x = np.zeros(len(coords))
    y = x.copy()
    dx = x.copy()
    dy = x.copy()
    for node in system.pot.G.nodes():
        xyz = system.node2xyz(node)
        i = system.pot.indices[node]
        x[i] = xyz[0]
        y[i] = xyz[1]
        
        theta = coords[i]
        dx[i] = np.cos(theta)
        dy[i] = np.sin(theta)
    
    if True:
        pot = system.get_potential()
        energies = pot.get_spin_energies(coords)
    else:
        energies = None
    
#    plt.quiver(x, y, dx, dy, energies, pivot="middle", cmap="YlOrRd_r")
    if energies is not None:
        plt.quiver(x, y, dx, dy, energies, pivot="middle", cmap="cool_r")
    else:
        plt.quiver(x, y, dx, dy, pivot="middle")
    d = 0.5
    plt.xlim([-d, dim[0]-d])
    plt.ylim([-d, dim[1]-d])
    plt.tight_layout()
    plt.xticks([])
    plt.yticks([])
    plt.gcf().set_facecolor('white')
    plt.box(on=False)

    plt.savefig(fname + ".eps", format="eps")
    

if __name__ == "__main__":
    main()