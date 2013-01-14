import ConfigParser
from os.path import expanduser

config = ConfigParser.ConfigParser()
config.add_section("gui")
config.set("gui", "use_pymol", value=False)

print "Open option file",expanduser("~")+'/.pygmin/pygmin.ini'
config.read(expanduser("~")+'/.pygmin/pygmin.ini')
