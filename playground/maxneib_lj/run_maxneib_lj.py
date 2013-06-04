from pele.potentials.maxneib_lj import MaxNeibsLJ, MaxNeibsLJSystem

def run_gui(system, db=None):
    import pele.gui.run as gr
    gr.run_gui(system, db=db)


if __name__ == "__main__":
    natoms = 20
    max_neibs=3.5
    rneib = 1.7
    epsneibs = 5.
    system = MaxNeibsLJSystem(natoms, max_neibs=max_neibs, rneib=rneib, epsneibs=epsneibs)

    dbname = "lj_N%d_n%.1f_rneib%.2f_epsn%.1f.db" %(natoms, max_neibs, rneib, epsneibs)
    print dbname

    
    run_gui(system, db=dbname)
