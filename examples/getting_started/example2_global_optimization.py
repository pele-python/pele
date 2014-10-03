"""
Global optimization
-------------------

We primarily use the global optimization method :class:`.Basinhopping` because it has been
shown to be very effective on the very high dimensional, smooth landscapes that we work with.

Basinhopping is iterative with each cycle composed of the following
features

1) random perturbation of the coordinates

2) local minimization

3) accept or reject the new coordinates based on the minimized function
   value

"""
p=1