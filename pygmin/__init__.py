import logging as _pygmin_logging
import sys as _pygmin_sys 

logger = _pygmin_logging.getLogger("pygmin")
global_handler = _pygmin_logging.StreamHandler(_pygmin_sys.stdout)
logger.addHandler(global_handler)
logger.setLevel(_pygmin_logging.DEBUG)

