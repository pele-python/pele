import ConfigParser
from os.path import expanduser

config = ConfigParser.ConfigParser()
config.add_section("gui")
config.set("gui", "use_pymol", value="False")

config.add_section("exec")
config.set("exec", "OPTIM", value="OPTIM", )
config.set("exec", "GMIN", value="GMIN")
config.set("exec", "AMBOPTIM", value="AMBOPTIM")

print "Open option file",expanduser("~")+'/.pele/pele.ini'
config.read(expanduser("~")+'/.pele/pele.ini')
