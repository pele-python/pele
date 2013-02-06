"""
tools to run the double ended connect in a separte process and 
make sure the the minima and transition states found are 
incorporated back into the master database
"""

import multiprocessing as mp
import threading as th
import sys
import time
from PyQt4 import QtCore, QtGui
import numpy as np

class UnboundMinimum(object):
    def __init__(self, minimum):
        self._id = minimum._id
        self.energy = minimum.energy
        self.coords = minimum.coords.copy()

class UnboundTransitionState(object):
    def __init__(self, ts):
        self._id = ts._id
        self.energy = ts.energy
        self.coords = ts.coords.copy()
        self.eigenvec = ts.eigenvec
        self.eigenval = ts.eigenval
        self._minimum1_id = ts._minimum1_id
        self._minimum2_id = ts._minimum2_id

class OutLog(object):
    """for redirecting stdout or stderr
    (edit, out=None, color=None) -> can write stdout, stderr to a
    QTextEdit.
    edit = QTextEdit
    out = alternate stream ( can be the original sys.stdout )
    color = alternate color (i.e. color stderr a different color)
    
    from http://www.riverbankcomputing.com/pipermail/pyqt/2009-February/022025.html
    """
    def __init__(self, conn):
        self.conn = conn

    def write(self, m):
        sys.stderr.write(":sending message:"+ m)
        self.conn.send(("stdout", m))
    
    def flush(self):
        pass


class DECProcess(mp.Process):
    """This object will run in a separate process and will actually do the connect run
    
    it will
    """
    def __init__(self, comm, system, min1, min2, pipe_stdout=True):
        mp.Process.__init__(self)
        #QtCore.QThread.__init__(self)
        self.comm = comm
        self.system = system
        self.min1, self.min2 = min1, min2
        self.pipe_stdout = pipe_stdout
    
    def clean_up(self):
        "send the lists of transition states and minima back to the parent process"
        minima = [UnboundMinimum(m) for m in self.db.minima()]
        tslist = [UnboundTransitionState(ts) for ts in self.db.transition_states()]
        self.comm.send(("new coords", minima, tslist))
    
    def do_double_ended_connect(self):
        db = self.system.create_database()
        self.db = db
        
        # min1 and min2 are associated with the old database, we need to create
        # the minima again using the new database
        m1 = db.addMinimum(self.min1.energy, self.min1.coords)
        m2 = db.addMinimum(self.min2.energy, self.min2.coords)
        
        connect = self.system.get_double_ended_connect(m1, m2, db,
                                                       fresh_connect=True)
        connect.connect()
    
    def run(self):
        if self.pipe_stdout:
#            print >> sys.stderr, "stderr"
            sys.stdout = OutLog(self.comm)
#            print >> sys.stderr, "stderr2"
        self.do_double_ended_connect()
        self.clean_up()

class Runner(QtCore.QObject):
    """this class will control spawning parallel jobs in the GUI"""
    def __init__(self, ProcessClass, *args, **kwargs):
        self.process_args = args
        self.process_kwargs = kwargs
        self.ProcessClass = ProcessClass
        self.mysubprocess = None
    
    def poll(self):
        print "polling"
        if not self.mysubprocess.is_alive():
            self.refresh_timer.stop()
            return
        if not self.parent_conn.poll():
            return
        
        message = self.parent_conn.recv()
        self.process_message(message)

    def start(self):
        if self.mysubprocess:
            if self.mysubprocess.is_alive():
                return
        parent_conn, child_conn = mp.Pipe()
        self.conn = parent_conn
        self.parent_conn = parent_conn
        
        self.refresh_timer = QtCore.QTimer()
        self.refresh_timer.timeout.connect(self.poll)
        self.refresh_timer.start(0.)

        self.mysubprocess = self.ProcessClass(child_conn, *self.process_args, **self.process_kwargs)
        self.mysubprocess.start()
    

    def process_message(self, message):
        raise NotImplementedError

class DECRunnerNEW(Runner):
    """this is the object that the gui creates to spawn a double ended connect run
    
    This will spawn a new process and deal with the communication
    """
    def __init__(self, system, database, min1, min2):
        super(DECRunner, self).__init__(DECProcess, system, min1, min2)
        self.database = database

    def add_minima_transition_states(self, new_minima, new_ts):
        print "processing new minima and ts"
        old2new = dict()
        for m in new_minima:
            mnew = self.database.addMinimum(m.energy, m.coords)
            old2new[m._id] = mnew
        
        for ts in new_ts:
            m1id = ts._minimum1_id
            m2id = ts._minimum2_id
            m1new = old2new[m1id]
            m2new = old2new[m2id]
            self.database.addTransitionState(ts.energy, ts.coords, m1new, 
                                             m2new, eigenval=ts.eigenval, 
                                             eigenvec=ts.eigenvec)
        nmin = len(new_minima)
        nts = len(new_ts)
        print "finished connect run: adding", nmin, "new minima, and", nts, "new transition states to database"

    
    def process_message(self, message):
        new_minima, new_ts = message
        self.add_minima_transition_states(new_minima, new_ts)


class DECRunner(QtCore.QObject):
    """this is the object that the gui creates to spawn a double ended connect run
    
    This will spawn a new process and deal with the communication
    """
    def __init__(self, system, database, min1, min2, outstream=None):
        QtCore.QObject.__init__(self)
        self.system = system
        self.database = database
        self.min1, self.min2 = min1, min2
        
        self.outstream = outstream
        
        #child_conn = self
        self.decprocess = None
#        self.lock = th.Lock()

    def poll(self):
#        if not self.decprocess.is_alive():
#            self.refresh_timer.stop()
#            return
        if not self.parent_conn.poll():
            return
        
        message = self.parent_conn.recv()
        self.process_message(message)
    
    def start(self):
        if(self.decprocess):
            if(self.decprocess.is_alive()):
                return
        parent_conn, child_conn = mp.Pipe()
        self.conn = parent_conn
        self.parent_conn = parent_conn
        
        self.decprocess = DECProcess(child_conn, self.system, self.min1, self.min2,
                                     pipe_stdout=(self.outstream is not None))
        self.decprocess.start()
    
#        self.poll_thread = PollThread(self, parent_conn)
#        self.poll_thread.start()
        self.refresh_timer = QtCore.QTimer()
        self.refresh_timer.timeout.connect(self.poll)
        self.refresh_timer.start(0.)


    def add_minima_transition_states(self, new_minima, new_ts):
        print "processing new minima and ts"
        old2new = dict()
        for m in new_minima:
            mnew = self.database.addMinimum(m.energy, m.coords)
            old2new[m._id] = mnew
        
        for ts in new_ts:
            m1id = ts._minimum1_id
            m2id = ts._minimum2_id
            m1new = old2new[m1id]
            m2new = old2new[m2id]
            self.database.addTransitionState(ts.energy, ts.coords, m1new, 
                                             m2new, eigenval=ts.eigenval, 
                                             eigenvec=ts.eigenvec)
        nmin = len(new_minima)
        nts = len(new_ts)
        print "finished connect run: adding", nmin, "new minima, and", nts, "new transition states to database"

    
    def process_message(self, message):
        sys.stderr.write("recieved message" + message[0] + "\n")
        if message[0] == "new coords":
            new_minima, new_ts = message[1:]
            self.add_minima_transition_states(new_minima, new_ts)
        elif message[0] == "stdout":
#            print "recieved stdout message"
#            print >> sys.stderr,  "recieved stdout message (stderr)", message[1]
            sys.stderr.write(":recieving message:" + message[1])
            self.outstream.write(message[1])
        
        