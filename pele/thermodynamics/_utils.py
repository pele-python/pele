"""
routines for computing thermodynamic information
"""
from __future__ import print_function
import multiprocessing as mp
import sys

from pele.thermodynamics._normalmodes import NormalModeError


class _ThermoWorker(mp.Process):  # pragma: no cover (coverage can't see it because it's in a separate process)
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

    def process_input(self):
        """get input from queue and process it
        
        return True if the queue is empty, raise any exception that occurs
        """
        # run until the input queue is empty
        if self.input_queue.empty():
            if self.verbose:
                print("worker ending")
            return True

        # get the next minima / ts to evaluate
        mts, mid, coords = self.input_queue.get()
        if mts == "ts":
            nnegative = 1
        # print "computing thermodynamics for ts", mid
        elif mts == "m":
            nnegative = 0
        else:
            raise Exception("mts must be 'm' or 'ts'")

        # do the computation
        invalid = False
        pgorder = self.system.get_pgorder(coords)
        try:
            fvib = self.system.get_log_product_normalmode_freq(coords, nnegative=nnegative)
        except NormalModeError as e:
            fvib = None
            invalid = True
            if mts == "m":
                sys.stdout.write(
                    "Problem computing normal modes for minimum with id {}. Setting m.invalid=True\n".format(mid))
            else:
                sys.stdout.write(
                    "Problem computing normal modes for transition state with id {}. Setting ts.invalid=True\n".format(
                        mid))
            sys.stdout.write(str(e) + "\n")
        if self.verbose:
            print("finished computing thermodynamic info for minimum", mid, pgorder, fvib)

        self.output_queue.put((mts, mid, fvib, pgorder, invalid))


    def run(self):
        while True:
            try:
                ret = self.process_input()
                if ret:
                    return
            except NormalModeError:
                # ignore these, they have already been dealt with
                continue
            except Exception as e:
                self.output_queue.put(e)


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
                 recalculate=False, commit_interval=100):
        self.system = system
        self.database = database
        self.verbose = verbose
        self.recalculate = recalculate
        self.commit_interval = commit_interval

        # initialize workers
        self.workers = []
        self.send_queue = mp.Queue()
        self.done_queue = mp.Queue()
        for _ in range(npar):
            worker = _ThermoWorker(self.send_queue, self.done_queue, system, verbose=self.verbose)
            worker.daemon = True
            self.workers.append(worker)

    def _populate_queue(self):
        """load the jobs into the queue
        """
        self.njobs = 0
        nmin = 0
        nts = 0
        for m in self.database.minima():
            if self.recalculate or (m.pgorder is None or m.fvib is None):
                self.njobs += 1
                self.send_queue.put(("m", m.id(), m.coords))
                nmin += 1

        for ts in self.database.transition_states():
            if self.recalculate or (ts.pgorder is None or ts.fvib is None):
                self.njobs += 1
                self.send_queue.put(("ts", ts.id(), ts.coords))
                nts += 1
        if self.verbose:
            print("computing thermodynamic info for {} minima and {} transition states".format(nmin, nts))

    def _process_return_value(self, ret):
        # if the a worker throws an unexpected exception, kill the workers and raise it
        if isinstance(ret, BaseException):
            self._kill_workers()
            raise ret
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

        # def _all_dead(self):
        #        for worker in self.workers:
        #            if worker.is_alive():
        #                return False
        #        return True

    def _get_results(self):
        """receive the results from the return queue
        """
        i = 0
        from queue import Empty
        while i < self.njobs:
            try:
                ret = self.done_queue.get()
            except Empty as e:
                sys.stderr.write("the queue is empty when it shouldn't be\n")
                self._kill_workers()
                raise e
            self._process_return_value(ret)
            if i % self.commit_interval == 0:
                self.database.session.commit()
                if self.verbose:
                    print("committing changes to the database")
            i += 1

    def finish(self):
        """kill the workers cleanly"""
        self.database.session.commit()
        if self.verbose:
            print("closing workers normally")
        for worker in self.workers:
            worker.join()
            worker.terminate()
            worker.join()

    def _kill_workers(self):
        """kill the workers without waiting for them to finish"""
        self.database.session.commit()
        if self.verbose:
            print("killing all workers")
        for worker in self.workers:
            worker.terminate()
            worker.join()


    def start(self):
        """start the computations
        
        this should be called after __init__
        """
        # populate the queue
        self._populate_queue()

        # start the workers
        for worker in self.workers:
            worker.start()

        # process the results as they come back
        self._get_results()

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
        print("computing fvib for minimum", m.id(), m.energy)
        m.fvib = system.get_log_product_normalmode_freq(m.coords)
    if commit:
        database.session.commit()
    return changed


def get_thermodynamic_information(system, database, nproc=4, recalculate=False, verbose=False):
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
        computer = GetThermodynamicInfoParallel(system, database, npar=nproc,
                                              recalculate=recalculate, verbose=verbose)
        computer.start()
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

def test():  # pragma: no cover
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

    # get_thermodynamic_information(system, db)

    print("getting thermodynamic info", db.number_of_minima())
    get_thermodynamic_information(system, db, nproc=4)

    for m in db.minima():
        print(m.id(), m.pgorder, m.fvib)

    print("\nnow transition states")
    for ts in db.transition_states():
        print(ts.id(), ts.pgorder, ts.fvib)


if __name__ == "__main__":
    test()

