"""
tools to run the double ended connect in a separte process and 
make sure the the minima and transition states found are 
incorporated back into the master database
"""

import multiprocessing as mp
import sys
import signal
import logging
import numpy as np

from PyQt4 import QtCore, QtGui

from pele.utils.events import Signal

class UnboundMinimum(object):
    def __init__(self, minimum):
        self._id = minimum._id
        self.energy = minimum.energy
        self.coords = minimum.coords.copy()

    def id(self):
        return self._id


class UnboundTransitionState(object):
    def __init__(self, ts):
        self._id = ts._id
        self.energy = ts.energy
        self.coords = ts.coords.copy()
        self.eigenvec = ts.eigenvec
        self.eigenval = ts.eigenval
        self._minimum1_id = ts._minimum1_id
        self._minimum2_id = ts._minimum2_id

    def id(self):
        return self._id


class OutLog(object):
    """for redirecting stdout or stderr
    
    everytime something is written to this object, it is sent through
    the pipe `conn`.  
    
    from http://www.riverbankcomputing.com/pipermail/pyqt/2009-February/022025.html
    """
    def __init__(self, conn):
        self.conn = conn
        self.message = ""

    def write(self, m):
        if len(m) > 0:
            self.conn.send(("stdout", m))
        return
##        sys.stderr.write(":sending message:"+ m)
#        self.message += m
##        if len(self.message) > 100:
##            self.flush()
###            self.conn.send(("stdout", m))
##        if len(self.mes)
#        if self.message[-1] == "\n":
#            self.flush()
    
    def flush(self):
#        self.conn.send(("stdout", self.message))
#        self.message = ""
        pass


class DECProcess(mp.Process):
    """This object will run in a separate process and will actually do the connect run
    
    when the run is finished the minima and transition states found will be sent
    back through the pipe as UnboundMinimum and UnboundTransitionState objects
    
    Parameters
    ----------
    comm : pipe
        child end of a mp.Pipe()
    system : 
    min1, min2 : 
        the minima to try to connect
    pipe_stdout : bool
        if true log messages will be sent back through the pipe
    return_smoothed_path : bool
        if the run ends successfully the smoothed path will be sent
        back through the pipe
    
    """
    def __init__(self, comm, system, min1, min2, pipe_stdout=True,
                  return_smoothed_path=True):
        mp.Process.__init__(self)
        #QtCore.QThread.__init__(self)
        self.comm = comm
        self.system = system
        self.min1, self.min2 = min1, min2
        self.pipe_stdout = pipe_stdout
        self.return_smoothed_path = return_smoothed_path
        
        self.started = False
        self.finished = False

    def get_smoothed_path(self):
        mints, S, energies = self.connect.returnPath()
        clist = [m.coords for m in mints]
        smoothpath = self.system.smooth_path(clist)
        return smoothpath, S, energies
    
    def test_success(self):
        return self.connect.graph.areConnected(self.m1local, self.m2local)
    
    def clean_up(self):
        """send the lists of transition states and minima back to the parent process"""
        minima = [UnboundMinimum(m) for m in self.db.minima()]
        tslist = [UnboundTransitionState(ts) for ts in self.db.transition_states()]
        self.comm.send(("new coords", minima, tslist))
        
        # return the success status
        success = self.test_success()
        self.comm.send(("success", success))
        
        if success:
            # return the smoothed path, or None if not successful
            pathdata = self.get_smoothed_path()
            self.comm.send(("smoothed path", pathdata))
        
        # send signal we're done here
        self.finished = True
        self.comm.send(("finished",))
    
    def terminate_early(self, *args, **kwargs):
        sys.stderr.write("caught signal, cleaning up and exiting\n")
        if self.started and not self.finished:
            sys.stderr.write("starting clean up\n")
            self.clean_up()
            sys.stderr.write("finished clean up\n")
        sys.stderr.write("exiting\n")
        sys.exit(0)
    
    def do_double_ended_connect(self):
        db = self.system.create_database()
        self.db = db
        
        # min1 and min2 are associated with the old database, we need to create
        # the minima again using the new database
        self.m1local = db.addMinimum(self.min1.energy, self.min1.coords)
        self.m2local = db.addMinimum(self.min2.energy, self.min2.coords)
        
        self.started = True
        self.connect = self.system.get_double_ended_connect(self.m1local, self.m2local, db,
                                                       fresh_connect=True)
        self.connect.connect()
    
    def run(self):
        signal.signal(signal.SIGTERM, self.terminate_early)
        signal.signal(signal.SIGINT, self.terminate_early)
        if self.pipe_stdout:
#            print >> sys.stderr, "stderr"
            self.mylog = OutLog(self.comm)
            sys.stdout = self.mylog
            logger = logging.getLogger("pele")
            handles = logger.handlers
            for h in handles:
#                print >> sys.stderr, "removing handler", h 
                logger.removeHandler(h)
            sh = logging.StreamHandler(self.mylog)
            logger.addHandler(sh)
#            import pele
#            logger.removeHandler(pele.h)
#            print >> sys.stderr, "stderr2"

        self.do_double_ended_connect()
        self.clean_up()


class DECRunner(QtCore.QObject):
    """Spawn a double ended connect run in a child process

    This will spawn a new process and deal with the communication
    
    
    Parameters
    ----------
    system : 
    database : Database
        The minima and transition states found will be added to the 
        database after the connect run is finished
    min1, min2 : Munimum objects
        the minima to try to connect
    outstream : an object with attribute `outstream.write(mystring)`
        the log messages from the connect run will be redirected here
    return_smoothed_path : bool
        if True the final smoothed path will be calculated
    
    Attributes
    ----------
    on_finished : Signal
        this signal will be called when the connect job is finished
    
    """
    def __init__(self, system, database, min1, min2, outstream=None,
                  return_smoothed_path=True, daemon=True):
        QtCore.QObject.__init__(self)
        self.system = system
        self.database = database
        self.min1, self.min2 = min1, min2
        self.return_smoothed_path = return_smoothed_path
        self.daemon = daemon
        
        self.outstream = outstream
        
        self.on_finished = Signal()
        
        self.decprocess = None
        
        self.newminima = set()
        self.newtransition_states = set()
        self.success = False
        self.killed_early = False
        self.is_running = False


    def poll(self):
        """this does the checking in the background to see if any messages have been passed"""
#        if not self.decprocess.is_alive():
#            self.refresh_timer.stop()
#            return
        if not self.parent_conn.poll():
            return
        
        message = self.parent_conn.recv()
        self.process_message(message)
    
    def start(self):
        """start the connect job"""
        if self.decprocess:
            if self.decprocess.is_alive():
                return
        parent_conn, child_conn = mp.Pipe()
        self.conn = parent_conn
        self.parent_conn = parent_conn
        
        self.decprocess = DECProcess(child_conn, self.system, self.min1, self.min2,
                                     pipe_stdout=(self.outstream is not None))
        self.decprocess.daemon = self.daemon
        self.decprocess.start()
    
#        self.poll_thread = PollThread(self, parent_conn)
#        self.poll_thread.start()
        self.refresh_timer = QtCore.QTimer()
        self.refresh_timer.timeout.connect(self.poll)
        self.refresh_timer.start(1.)
        self.is_running = True


    def add_minima_transition_states(self, new_minima, new_ts):
        """Add the minima and transition states found to the database
        
        convert the UnboundMinimum and UnboundTransitionStates to ones
        bound to self.database
        """
        print "processing new minima and ts"
        self.newminima = set()
        self.newtransition_states = set()
        old2new = dict()
        self.system.params.gui._sort_lists = False
        for m in new_minima:
            mnew = self.database.addMinimum(m.energy, m.coords)
            old2new[m._id] = mnew
            self.newminima.add(mnew)
        
        for ts in new_ts:
            m1id = ts._minimum1_id
            m2id = ts._minimum2_id
            m1new = old2new[m1id]
            m2new = old2new[m2id]
            tsnew = self.database.addTransitionState(ts.energy, ts.coords, m1new, 
                                             m2new, eigenval=ts.eigenval, 
                                             eigenvec=ts.eigenvec)
            self.newtransition_states.add(tsnew)
        nmin = len(new_minima)
        nts = len(new_ts)
        print "finished connect run: adding", nmin, "minima, and", nts, "transition states to database"
        self.system.params.gui._sort_lists = True

    def terminate_early(self):
        self.killed_early = True
        self.decprocess.terminate()
        print "finished terminating"
        self.is_running = False
#        self.decprocess.join()
#        print "done killing job"
#        self.on_finished()
    
    def finished(self):
        """the job is finished, do some clean up"""
        self.decprocess.join()
        self.decprocess.terminate()
        self.decprocess.join()
        self.refresh_timer.stop()
#        print "done killing job"
        self.on_finished()
        self.is_running = False
        
    
    def process_message(self, message):
        if message[0] == "stdout":
            self.outstream.write(message[1])
        elif message[0] == "new coords":
            new_minima, new_ts = message[1:]
            self.add_minima_transition_states(new_minima, new_ts)
        elif message[0] == "success":
            self.success = message[1]
        elif message[0] == "smoothed path":
            pathdata = message[1]
            self.smoothed_path, self.S, self.energies = pathdata
        elif message[0] == "finished":
            self.finished()
        
        
