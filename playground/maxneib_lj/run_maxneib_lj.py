from pygmin.potentials.maxneib_lj import MaxNeibsLJ, MaxNeibsLJSystem

def run_gui(system, db=None):
    import pygmin.gui.run as gr
    gr.run_gui(system)


if __name__ == "__main__":
    natoms = 20
    system = MaxNeibsLJSystem(natoms, max_neibs=3, rneib=1.7, epsneibs=5.)
    
    run_gui(system)# , db="maxneib_lj.db")
