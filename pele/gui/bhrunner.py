from __future__ import print_function
import time
import multiprocessing as mp

import numpy as np
from PyQt5 import QtCore

from pele.utils.events import Signal

class _BHProcess(mp.Process):
    """do basinhopping in a different process
    
    this class actually does the basinhopping run and the minima it finds to 
    it's parent through a communication pipe
    """
    def __init__(self, system, comm, nsteps=None):
        mp.Process.__init__(self)
        self.comm = comm
        self.system = system
        self.nsteps = nsteps
    
    def run(self):
        seed = int(time.time() * 100.) % 4294967295
        np.random.seed(seed)
        print(np.random.random(2))
        db = self.system.create_database()
        db.on_minimum_added.connect(self.insert)
        opt = self.system.get_basinhopping(database=db, outstream=None)
        if self.nsteps is None:
            try:
                self.nsteps = self.system.params.gui.basinhopping_nsteps
            except AttributeError or KeyError:
                self.nsteps = 100
        
        for i in range(self.nsteps):
            opt.run(1)
            if self._should_die():
                return
    
    def _should_die(self):
        """check if we have received a message telling us to die"""
        if self.comm.poll():
            message = self.comm.recv()
            if message == "kill":
                return True
            else:
                print("don't understand message", message, "ignoring")
        return False
       
    def insert(self, m):
        self.comm.send([m.energy,m.coords])
        
class PollThread(QtCore.QThread):
    def __init__(self, bhrunner, conn):
        QtCore.QThread.__init__(self)
        self.bhrunner = bhrunner
        self.conn = conn
    def run(self):
        while self.bhrunner.bhprocess.is_alive():
            if self.conn.poll():
                minimum = self.conn.recv()
                self.emit(QtCore.SIGNAL("Activated( PyQt_PyObject )"),minimum)                


class BHRunner(QtCore.QObject):
    """manage a single basinhopping run in a separate process
    
    This class spawns the basinhopping job in a separate process and receives
    the minima found through a pipe
    """
    def __init__(self, system, database, nsteps=None, on_finish=None):
        QtCore.QObject.__init__(self)
        self.system = system
        self.database = database
        self.daemon = True
        self.nsteps = nsteps
        
        #child_conn = self
        self.bhprocess = None
        
        self.on_finish = Signal()
        if on_finish is not None:
            self.on_finish.connect(on_finish)
    
    def is_alive(self):
        return self.bhprocess.is_alive()
    
    def poll(self):
        if not self.parent_conn.poll():
            if not self.bhprocess.is_alive():
                self.refresh_timer.stop()
                self.on_finish()
                return
            return
        
        minimum = self.parent_conn.recv()
        self.database.addMinimum(minimum[0],minimum[1])        

    def start(self):
        if self.bhprocess:
            if self.bhprocess.is_alive():
                return
        parent_conn, child_conn = mp.Pipe()
        
        self.bhprocess = _BHProcess(self.system, child_conn, nsteps=self.nsteps)
        self.bhprocess.daemon = self.daemon
        self.bhprocess.start()
        self.parent_conn = parent_conn
        self.refresh_timer = QtCore.QTimer()
        self.refresh_timer.timeout.connect(self.poll)
        self.refresh_timer.start(50.) # time in msec
    
    def kill(self):
        """kill the job that is running"""
        if self.is_alive():
            self.parent_conn.send("kill")
            self.bhprocess.join()

class BHManager(object):
    def __init__(self, system, database, on_number_alive_changed=None):
        self.system = system
        self.database = database
        self.workers = []
        self._old_nalive = -1
        
        self.on_number_alive_changed = Signal()
        if on_number_alive_changed is not None:
            self.on_number_alive_changed.connect(on_number_alive_changed)

        self.refresh_timer = QtCore.QTimer()
        self.refresh_timer.timeout.connect(self._check_number)

        
    def _remove_dead(self):
        self.workers = [w for w in self.workers if w.is_alive()]
        if self._old_nalive != len(self.workers):
            self._old_nalive = len(self.workers)
            self.on_number_alive_changed(len(self.workers))
    
    def _check_number(self):
        self._remove_dead()
        if len(self.workers) == 0:
            self.refresh_timer.stop()
    
    def start_worker(self, nsteps=None):
        self._remove_dead()
        worker = BHRunner(self.system, self.database, nsteps=nsteps)
        worker.on_finish.connect(self._remove_dead)
        worker.start()
        self.workers.append(worker)
    
        if not self.refresh_timer.isActive():
            self.refresh_timer.start(1000.) # time in msec

    
    def kill_all_workers(self):
        for worker in self.workers:
            worker.kill()
    
    def number_of_workers(self):
        self._remove_dead()
        return len(self.workers)
        
def found(m):
    print(("New Minimum",m[0]))
    
if __name__ == '__main__': 
    runner = BHRunner(onMinimumAdded=found)
    runner.start()
    #bhp.start()
    #while(True):
    #    pass# bhp.is_alive():
        #print "is alive"
    #bhp.join()

