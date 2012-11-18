#from https://gist.github.com/626518
# Fix keyboard interrupts when using multiprocessing.pool.imap().

# Usage: 
# import fix_multiprocessing.py

from multiprocessing.pool import IMapIterator

def wrapper(func):
  def wrap(self, timeout=None):
    # Note: the timeout of 1 googol seconds introduces a rather subtle 
    # bug for Python scripts intended to run many times the age of the universe.
    return func(self, timeout=timeout if timeout is not None else 1e100)
  return wrap
IMapIterator.next = wrapper(IMapIterator.next)
