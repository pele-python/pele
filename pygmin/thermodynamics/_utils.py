
     
    
def get_thermodynamic_information_minimum(system, database, minimum, commit=True):
    m = minimum
    changed = False
    if m.pgorder is None:
        changed = True
        m.pgorder = system.get_pgorder(m.coords)
    if m.fvib is None:
        changed = True
        print "computing fvib for minima", m._id, m.energy
        m.fvib = system.get_log_product_normalmode_freq(m.coords)
    if commit:
        database.session.commit()
    return changed


def get_thermodynamic_information(system, database):
    """
    compute thermodynamic information for all minima in a database
    
    Parameters
    ----------
    system : pygmin System class
    databse : a Database object
    
    Notes
    -----
    The information that is computed is the point group order (m.pgorder) and the
    log product of the squared normal mode frequencies (m.fvib).
    """
    changed = False
    try:
        for m in database.minima():
            c = get_thermodynamic_information_minimum(system, database, m, commit=False)
            if c: changed = True
    except KeyboardInterrupt:
        if changed:
            database.session.commit()
        raise

    if changed:    
        database.session.commit()
    
