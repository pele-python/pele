from optparse import OptionParser
from pygmin.storage.database import Database
import sys

class Application(object):   
    options = None
    compareStructures = None
    database = None
    parser = OptionParser()
    option_groups = {}
    
    def __init__(self):
        pass
    
    def create_potential(self):
        raise Exception("system does not implement create_potential")
    
    def add_options(self):
        pass
    
    def add_option(self, *args, **kwargs):
        group = None
        if kwargs.has_key("group"):
            del kwargs["group"]
        grp = self.parser
        if not group is None:
            if not self.option_groups.has_key(group):
                grp = self.parser.add_option_group(group)
                self.option_groups[group] = grp
        grp.add_option(*args, **kwargs)
    
    def run(self):
        raise Exception("run not implemented for application")
    
    def execute(self):
        print "The program was started with the following command line argument"
        print " ".join(sys.argv)
        self.add_options()        
        (self.options, self.args) = self.parser.parse_args()
        print "Starting the run with the following options"
        print "-------------------------------------------"
        for key, value in self.options.__dict__.iteritems():
            print str(key) + " = " + str(value)
        print "-------------------------------------------"
        self.database = Database(db=self.options.database, accuracy=self.accuracy, compareMinima=self.compareStructures)
        self.run()

                
    