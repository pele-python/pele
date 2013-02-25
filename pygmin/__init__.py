import logging
import sys

logger = logging.getLogger("pygmin")
h = logging.StreamHandler(sys.stdout)
logger.addHandler(h)
logger.setLevel(logging.DEBUG)

