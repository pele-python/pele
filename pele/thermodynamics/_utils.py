"""
routines for computing thermodynamic information
"""

import multiprocessing as mp
from pele.thermodynamics._normalmodes import NormalModeError

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
    def __init__(self, input_queue, output_queue, system, verbose=False, **kwargs):
        mp.Process.__init__(self, **kwargs)
        self.input_queue = input_queue
        self.output_queue = output_queue
        self.system = system
        self.verbose = verbose


    def run(self):
        while True:
            # run until the input queue is empty
            if self.input_queue.empty():
                if self.verbose:
                    print "woker ending"
                return
            
            # get the next minima / ts to evaluate
            mts, mid, coords = self.input_queue.get()
            if mts == "ts":
                nnegative=1
#                print "computing thermodynamics for ts", mid
            elif mts == "m":
                nnegative=0
            else:
                raise Exception("mts must be 'm' or 'ts'")
            
            # do the computation
            invalid = False
            pgorder = self.system.get_pgorder(coords)
            try:
                fvib = self.system.get_log_product_normalmode_freq(coords, nnegative=nnegative)
            except NormalModeError, e:
                fvib = None
                invalid = True
                if mts == "m":
                    print "problem computing normal modes for minimum with id", mid
                else:  
                    print "problem computing normal modes for transition state with id", mid
                print str(e)
            if self.verbose:
                print "finished computing thermodynamic info for minimum", mid, pgorder, fvib
            
            self.output_queue.put((mts, mid, fvib, pgorder, invalid))


class GetThermodynamicInfoParallel(object):
    """
    a class to compute thermodynamic information in parallel
    
    Parameters
    ----------
    system : pele system object
    database : pele database
        thermodynamic information will be calculated for all minima
        in the database that don't already have the data
    verbose : bool
        specify verbosity
    only_minima : bool
        if True the transition state free energy will not be computed
    """
    def __init__(self, system, database, npar=4, verbose=False, only_minima=False, 
                  recalculate=False):
        self.system = system
        self.database = database
        self.verbose = verbose
        self.recalculate = recalculate

        # initialize workers
        self.workers = []
        self.send_queue = mp.Queue()
        self.done_queue = mp.Queue()
        for i in range(npar):
            worker = _ThermoWorker(self.send_queue, self.done_queue, system, verbose=self.verbose)
            worker.daemon = True
            self.workers.append(worker)
    
    def _populate_queue(self):
        """load the jobs into the queue
        """
        self.njobs = 0
        for m in self.database.minima():
            if self.recalculate or (m.pgorder is None or m.fvib is None):
                self.njobs += 1
                self.send_queue.put(("m", m._id, m.coords))
        
        for ts in self.database.transition_states():
            if self.recalculate or (ts.pgorder is None or ts.fvib is None):
                self.njobs += 1
                self.send_queue.put(("ts", ts._id, ts.coords))
    
    def _process_return_value(self, ret):
        mts, mid, fvib, pgorder, invalid = ret
        if mts == "m":
            m = self.database.getMinimum(mid)
            m.fvib = fvib
            m.pgorder = pgorder
            if invalid:
                m.invalid = True
        elif mts == "ts":
            ts = self.database.getTransitionStateFromID(mid)
            ts.fvib = fvib
            ts.pgorder = pgorder
            if invalid:
                ts.invalid = True
        else:
            raise Exception("mts must be 'm' or 'ts'")
        self.database.session.commit()

      
    def _get_results(self):
        """receive the results from the return queue
        """
        for i in xrange(self.njobs):
            ret = self.done_queue.get()
            self._process_return_value(ret)
    
    def finish(self):
        if self.verbose:
            print "closing workers"
        for worker in self.workers:
            worker.join()
            worker.terminate()
            worker.join()



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
        self.finish()
        
    
    
    

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


def get_thermodynamic_information(system, database, nproc=4, recalculate=False):
    """
    compute thermodynamic information for all minima and transition states in a database
    
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
        worker = GetThermodynamicInfoParallel(system, database, npar=nproc, 
                                              recalculate=recalculate)
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
    from pele.landscape import ConnectManager
    system = LJCluster(13)
    
    db = system.create_database()
    bh = system.get_basinhopping(db, outstream=None)
    bh.run(200)
    
    manager = ConnectManager(db)
    for i in range(3):
        min1, min2 = manager.get_connect_job()
        connect = system.get_double_ended_connect(min1, min2, db)
        connect.connect()
    
#    get_thermodynamic_information(system, db)
    
    print "getting thermodynamic info", db.number_of_minima()
    get_thermodynamic_information(system, db, nproc=4)
    
    for m in db.minima():
        print m._id, m.pgorder, m.fvib
    
    print "\nnow transition states"
    for ts in db.transition_states():
        print ts._id, ts.pgorder, ts.fvib

if __name__ == "__main__":
    test()
