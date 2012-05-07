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
import threading as mp
#import multiprocessing as mp
import threading as th
import time

class BHProcess(mp.Thread): #Process):
    def __init__(self, comm):
        mp.Thread.__init__(self)
        self.comm = comm
    
    def run(self):
        import numpy as np
        import potentials.lj as lj
        import basinhopping as bh
        import take_step.random_displacement as random_displacement
        natoms = 13 

        # random initial coordinates
        coords = np.random.random(3 * natoms)
        potential = lj.LJ()#1.0, 1.0, None)
        step = random_displacement.takeStep(stepsize=0.5)
        opt = bh.BasinHopping(coords,potential,
                          temperature=5., takeStep=step.takeStep, storage=self.insert)
        print "start bh"
        #while(True):
        opt.run(2000)
        print "done with bh"
        
    def insert(self, E, coords):
        self.comm.send([E,coords])
        
class BHRunner():
    def __init__(self, onNewMinimumFound=None):
        self.storage = savenlowest.SaveN(50, onNewMinimumFound=onNewMinimumFound)
        
        #self.poll_thread = th.Thread(target=self.poll)
        
        #self.conn, child_conn = mp.Pipe()
        child_conn = self
        self.bhprocess = BHProcess(child_conn)
    
    def send(self, minimum):
        self.minimum_found(minimum[0], minimum[1])
    
    def pause(self):
        pass
    
    def contiue(self):
        pass
    
    def start(self):
        #self.poll_thread.start()
        self.bhprocess.start()
        pass
    
    def poll(self):    
        while(self.bhprocess.is_alive()):
            if(self.conn.poll()):
                minimum = self.conn.recv()
                self.minimum_found(minimum[0], minimum[1])
    
    def minimum_found(self,E,coords):
        self.storage.insert(E,coords)
    
def found(E,coords):
    print("New Minimum",E)
    
if __name__ == '__main__': 
    runner = BHRunner(onNewMinimumFound=found)
    runner.start()
    #bhp.start()
    #while(True):
    #    pass# bhp.is_alive():
        #print "is alive"
    #bhp.join()
