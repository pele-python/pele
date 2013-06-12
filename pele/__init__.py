import logging as _pele_logging
import sys as _pele_sys 

logger = _pele_logging.getLogger("pele")
global_handler = _pele_logging.StreamHandler(_pele_sys.stdout)
logger.addHandler(global_handler)
logger.setLevel(_pele_logging.DEBUG)

