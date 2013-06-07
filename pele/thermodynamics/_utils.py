"""
routines for computing thermodynamic information
"""

import multiprocessing as mp
from collections import deque

class _ThermoWorker(mp.Process):
    """worker to calculate the thermodynamic data in a separate process
    
    Parameters
    ----------
    input_queue : mp.Queue object
        the worker will receive jobs on this queue
    output_queue : mp.Queue object
        the worker will return results on this queue
    system : pele system object
    """
    def __init__(self, input_queue, output_queue, system, **kwargs):
        mp.Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.output_queue = output_queue
        self.system = system


    def run(self):
        while True:
            if self.input_queue.empty():
                return
            mid, coords = self.input_queue.get()
            pgorder = self.system.get_pgorder(coords)
            fvib = self.system.get_log_product_normalmode_freq(coords)
            
            self.output_queue.put((mid, fvib, pgorder))


class GetThermodynamicInfoParallel(object):
    """
    a class to compute thermodynamic information in parallel
    
    Parameters
    ----------
    system : pele system object
    database : pele database
        thermodynamic information will be calculated for all minima
        in the database that don't already have the data
    """
    def __init__(self, system, database, npar=4):
        self.system = system
        self.database = database
        # initialize workers
        self.workers = []
        self.send_queue = mp.Queue()
        self.done_queue = mp.Queue()
        for i in range(npar):
            worker = _ThermoWorker(self.send_queue, self.done_queue, system)
            worker.daemon = True
            self.workers.append(worker)
    
    def _populate_queue(self):
        """load the jobs into the queue
        """
        self.njobs = 0
        for m in self.database.minima():
            if m.pgorder is None or m.fvib is None:
                self.njobs += 1
                self.send_queue.put((m._id, m.coords))
        
    def _get_results(self):
        """receive the results from the return queue
        """
        for i in xrange(self.njobs):
            mid, fvib, pgorder = self.done_queue.get()
            m = self.database.getMinimum(mid)
            m.fvib = fvib
            m.pgorder = pgorder
            self.database.session.commit()
#            print "got result", mid
            
    def start(self):
        # populate the queue
        self._populate_queue()

        # start the workers
        for worker in self.workers:
            worker.start()
        
        # process the results as they come back
        try:    
            self._get_results()
        except:
            # kill all the child processes
            self.database.session.commit()
            for worker in self.workers:
                worker.terminate()
                worker.join()
            raise
        
        # kill the workers cleanly
        for worker in self.workers:
                worker.join()
                worker.terminate()
                worker.join()

        
    
    
    

def get_thermodynamic_information_minimum(system, database, minimum, commit=True):
    m = minimum
    changed = False
    if m.pgorder is None:
        changed = True
        m.pgorder = system.get_pgorder(m.coords)
    if m.fvib is None:
        changed = True
        print "computing fvib for minimum", m._id, m.energy
        m.fvib = system.get_log_product_normalmode_freq(m.coords)
    if commit:
        database.session.commit()
    return changed


def get_thermodynamic_information(system, database, nproc=None):
    """
    compute thermodynamic information for all minima in a database
    
    Parameters
    ----------
    system : pele System class
    database : a Database object
    nproc : number of processors to use
    
    Notes
    -----
    The information that is computed is the point group order (m.pgorder) and the
    log product of the squared normal mode frequencies (m.fvib).
    """
    if nproc is not None:
        worker = GetThermodynamicInfoParallel(system, database, npar=nproc)
        worker.start()
        return
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

#
# only testing stuff below here
#

def test():
    from pele.systems import LJCluster
    system = LJCluster(13)
    
    db = system.create_database()
    bh = system.get_basinhopping(db, outstream=None)
    bh.run(200)
    
#    get_thermodynamic_information(system, db)
    
    print "getting thermodynamic info", db.number_of_minima()
    get_thermodynamic_information(system, db, nproc=4)
    
    for m in db.minima():
        print m._id, m.pgorder, m.fvib
    

if __name__ == "__main__":
    test()
