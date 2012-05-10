# Comment: multiprocess not running yet, still using threading module
#
#
# All calculations need to be done without blocking the
# gui. The basin hopping runner takes care for running 
# basin hopping simulations. It uses subproceses instead of
# threading since threading is only possile on GIL level
# and will not work properly if native code is involved
#
# The layout is the following: The BHRunner in master has
# the storage class. The subprocess communicates all found
# minima to the master via pipe which handles the insert and
# gui updates.

from storage import savenlowest
#import threading as mp
import multiprocessing as mp
import threading as th
import time
from PyQt4 import QtCore,QtGui

class BHProcess(mp.Process):
    def __init__(self, comm):
        mp.Process.__init__(self)
        #QtCore.QThread.__init__(self)
        self.comm = comm
    
    def run(self):
        import numpy as np
        import potentials.lj as lj
        import basinhopping as bh
        import take_step.random_displacement as random_displacement
        natoms = 114 

        # random initial coordinates
        coords = np.random.random(3 * natoms)
        potential = lj.LJ()#1.0, 1.0, None)
        step = random_displacement.takeStep(stepsize=0.5)
        opt = bh.BasinHopping(coords,potential,
                          temperature=1., takeStep=step.takeStep, storage=self.insert)
        print "start bh"
        #while(True):
        opt.run(500)
        print "done with bh"
        
    def insert(self, E, coords):
            self.comm.send([E,coords])
        
class PollThread(QtCore.QThread):
    def __init__(self, bhrunner, conn):
        QtCore.QThread.__init__(self)
        self.bhrunner = bhrunner
        self.conn = conn
    def run(self):
        while(self.bhrunner.bhprocess.is_alive()):
            #print "no data"
            if(self.conn.poll()):
                minimum = self.conn.recv()
                self.emit(QtCore.SIGNAL("Activated( PyQt_PyObject )"),minimum)                


class BHRunner(QtCore.QObject):
    def __init__(self, onMinimumAdded=None, onMinimumRemoved=None):
        QtCore.QObject.__init__(self)
        self.storage = savenlowest.SaveN(100, 
                         onMinimumAdded=onMinimumAdded,
                         onMinimumRemoved=onMinimumRemoved)
        
        
        #child_conn = self
        self.bhprocess = None
        self.lock = th.Lock()
        
    def send(self, minimum):
        self.minimum_found(minimum)
    
    def pause(self):
        pass
    
    def contiue(self):
        pass
    
    def start(self):
        if(self.bhprocess):
            if(self.bhprocess.is_alive()):
                return
        parent_conn, child_conn = mp.Pipe()
        
        self.bhprocess = BHProcess(child_conn)
        self.bhprocess.start()
        self.poll_thread = PollThread(self, parent_conn)
        self.connect(self.poll_thread,QtCore.SIGNAL("Activated ( PyQt_PyObject ) "), self.minimum_found)        
        self.poll_thread.start()  
    def Activated(self, str):
        
        print str          
        
    def minimum_found(self,minimum):
        self.lock.acquire()
        self.storage.insert(minimum[0],minimum[1])
        self.lock.release()
    
def found(m):
    print("New Minimum",m[0])
    
if __name__ == '__main__': 
    runner = BHRunner(onMinimumAdded=found)
    runner.start()
    #bhp.start()
    #while(True):
    #    pass# bhp.is_alive():
        #print "is alive"
    #bhp.join()
