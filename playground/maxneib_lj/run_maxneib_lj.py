from pygmin.potentials.maxneib_lj import MaxNeibsLJ, MaxNeibsLJSystem

def run_gui(system, db=None):
    import pygmin.gui.run as gr
    gr.run_gui(system)


if __name__ == "__main__":
    natoms = 10
    system = MaxNeibsLJSystem(natoms, max_neibs=6, rneib=1.7)
    
    run_gui(system)# , db="maxneib_lj.db")
