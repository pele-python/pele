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

#import threading as mp
import multiprocessing as mp
import threading as th
import time
from PyQt4 import QtCore,QtGui
import numpy as np

class BHProcess(mp.Process):
    def __init__(self, system, comm):
        mp.Process.__init__(self)
        #QtCore.QThread.__init__(self)
        self.comm = comm
        self.system = system
    
    def run(self):
        seed = int(time.time())#*100.)
        print seed
        np.random.seed(seed)
        print np.random.random(2)
        opt = self.system.createBasinHopping()
        opt.storage = self.insert
        #while(True):
        print 'bhrunner.py: number of BH steps set to 1'
        opt.run(1)
        
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
    def __init__(self, system, onMinimumAdded=None, onMinimumRemoved=None):
        QtCore.QObject.__init__(self)
        
        self.system = system
        
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
        
        self.bhprocess = BHProcess(self.system, child_conn)
        self.bhprocess.start()
        self.poll_thread = PollThread(self, parent_conn)
        self.connect(self.poll_thread,QtCore.SIGNAL("Activated ( PyQt_PyObject ) "), self.minimum_found)        
        self.poll_thread.start()  
        
    def minimum_found(self,minimum):
        self.lock.acquire()
        self.system.storage.addMinimum(minimum[0],minimum[1])
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
