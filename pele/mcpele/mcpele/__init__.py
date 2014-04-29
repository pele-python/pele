import logging as _mcpele_logging
import sys as _mcpele_sys 

logger = _mcpele_logging.getLogger("mcpele")
global_handler = _mcpele_logging.StreamHandler(_mcpele_sys.stdout)
logger.addHandler(global_handler)
logger.setLevel(_mcpele_logging.DEBUG)

