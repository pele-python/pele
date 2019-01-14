"""
.. currentmodule:: pele.concurrent

Parallel connect jobs (`pele.concurrent`)
===========================================
Concurrent programming concepts allow pele connect jobs to be run in parallel. We use the package
Pyro4 which makes this amazingly simple.

The server manages the central database and decides which connections to try.

.. autosummary::
   :toctree: generated/

    ConnectServer

An arbitrary number of workers can connect to the server and process the connect jobs.

.. autosummary::
   :toctree: generated/

    ConnectWorker


Usage
-----
see the example in the examples/connecting_in_parallel/ folder for more details.

start the server in one terminal::

    $ python start_server.py

The provides an uri (which is also stored in pyros.uri). Clients can connect
to this uri and access the connect manager. To start a worker::

    $ python start_worker.py

Run on cluster / with remote workers
------------------------------------
Start the server on a workstation (or node) which should be the master node.
To allow for incoming remote connections, a hostname must be specified. Then
specify the hostname to connect to in worker.py and submit as many jobs a needed.


"""
from __future__ import absolute_import

from ._connect_server import *

