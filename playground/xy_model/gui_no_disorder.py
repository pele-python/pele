from xy_model_system import XYModlelSystem

def run_gui_nodisorder(L=24):
    from pele.gui import run_gui
    dim=[L,L]
    system = XYModlelSystem(dim=dim, phi_disorder=0.)
    system.params.basinhopping.temperature=10.
    db = system.create_database("xy_%dx%d_nodisorder.sqlite" %(L,L))
    
    run_gui(system, db=db)


run_gui_nodisorder(L=12)