from __future__ import division
import numpy as np
import abc
import tensorflow as tf
from pele.potentials import BasePotential
from itertools import permutations, combinations
from pele.transition_states._zeroev import orthogonalize
from scipy.special import factorial
from pele.potentials import MeanFieldPSpinSpherical
from numba import jit
from pele.optimize._quench import lbfgs_cpp
from pele.optimize import Result

def dot(x, y):
    return tf.reduce_sum(tf.mul(x,y))

def rnorm(x):
    return tf.rsqrt(tf.reduce_sum(tf.square(x)))

def orthog(x, zerov):
    return tf.sub(x, tf.mul(dot(x, zerov), zerov))

def MeanFieldPSpinSphericalTF(interactions, nspins, p, dtype='float64'):
    if p == 3:
        return MeanField3SpinSphericalTF(interactions, nspins, dtype=dtype)
    elif p == 4:
        return MeanField4SpinSphericalTF(interactions, nspins, dtype=dtype)
    elif p == 5:
        return MeanField5SpinSphericalTF(interactions, nspins, dtype=dtype)
    else:
        raise Exception("BaseMeanFieldPSpinSphericalTF: p={} not implemented".format(p))

class BaseMeanFieldPSpinSphericalTF(BasePotential):
    """
    the potential has been hardcoded for p=3
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, interactions, nspins, p, dtype='float64'):
        self.dtype = 'float64'
        self.p = p
        self.nspins = nspins
        prf = np.power(np.sqrt(self.nspins), self.p-1) if self.p > 2 else 1.
        self.prf = tf.Variable(tf.constant(prf, dtype=self.dtype))
        self.sqrtN = np.sqrt(self.nspins)
        interactions = self._adaptInteractions(interactions)
        self.interactions = tf.Variable(tf.constant(interactions, dtype=self.dtype))
        self.x1 = tf.Variable(tf.zeros([self.nspins], dtype=self.dtype))
        self.x2 = tf.stop_gradient(self.x1) #this makes an efficient copy of the variable x1 without it contributing to the gradient
        self.session = tf.Session()
        init = tf.initialize_all_variables()
        self.session.run(init)

    def __exit__(self, exc_type, exc_value, traceback):
        self.session.close()

    def _adaptInteractions(self, interactions):
        """
        this function guarantees that the repeated indexes terms are zero
        :param nspins:
        :param p:
        :return:
        """
        new_interactions = np.zeros(interactions.shape)
        for comb in combinations(range(self.nspins), self.p):
                for perm in permutations(comb):
                    new_interactions[perm] = interactions[perm]
        return new_interactions

    @property
    def lossTensorPartial(self):
        """this gives a 2x2 matrix generator"""
        return tf.mul(tf.reshape(self.x1, [self.nspins]), tf.reshape(self.x2, [self.nspins,1]))

    @property
    @abc.abstractmethod
    def lossTensor(self):
        """this gives the full loss tensor generator"""

    @property
    def loss(self):
        """this mutliplies times the interaction and reduces sum reduces the tensor"""
        return -tf.div(tf.reduce_sum(tf.mul(self.interactions, self.lossTensor)), self.prf)

    @property
    def compute_gradient(self):
        grad = tf.gradients(self.loss, self.x1)[0]
        e = tf.reduce_sum(tf.mul(self.x1, grad))
        grad = orthog(grad, tf.mul(self.x1, rnorm(self.x1)))
        return e, grad

    @jit
    def _normalizeSpins(self, coords):
        coords /= (np.linalg.norm(coords)/self.sqrtN)

    def getEnergy(self, coords):
        self._normalizeSpins(coords)
        return self.session.run(self.loss, feed_dict={self.x1 : coords})

    def getEnergyGradient(self, coords):
        self._normalizeSpins(coords)
        e, grad = self.session.run(self.compute_gradient,
                                   feed_dict={self.x1 : coords})
        return e, grad

class MeanField3SpinSphericalTF(BaseMeanFieldPSpinSphericalTF):
    def __init__(self, interactions, nspins, dtype='float64'):
        super(MeanField3SpinSphericalTF, self).__init__(interactions, nspins, p=3, dtype=dtype)

    @property
    def lossTensor(self):
        """this gives the full loss tensor generator"""
        return tf.mul(self.lossTensorPartial, tf.reshape(self.x2, [self.nspins,1,1]))

class MeanField4SpinSphericalTF(BaseMeanFieldPSpinSphericalTF):
    def __init__(self, interactions, nspins, dtype='float64'):
        super(MeanField4SpinSphericalTF, self).__init__(interactions, nspins, p=4, dtype=dtype)

    @property
    def lossTensor(self):
        """this gives the full loss tensor generator"""
        return tf.mul(tf.mul(self.lossTensorPartial, tf.reshape(self.x2, [self.nspins,1,1])),
                      tf.reshape(self.x2, [self.nspins,1,1,1]))

class MeanField5SpinSphericalTF(BaseMeanFieldPSpinSphericalTF):
    def __init__(self, interactions, nspins, dtype='float64'):
        super(MeanField5SpinSphericalTF, self).__init__(interactions, nspins, p=5, dtype=dtype)

    @property
    def lossTensor(self):
        """this gives the full loss tensor generator"""
        return tf.mul(tf.mul(tf.mul(self.lossTensorPartial, tf.reshape(self.x2, [self.nspins,1,1])),
                             tf.reshape(self.x2, [self.nspins,1,1,1])),
                      tf.reshape(self.x2, [self.nspins,1,1,1,1]))

class GradientDescent(tf.train.GradientDescentOptimizer):
    def __init__(self, potential, learning_rate=tf.constant(0.1, tf.float64), use_locking=False, tol=1e-5):
        super(GradientDescent, self).__init__(learning_rate, use_locking=use_locking, name='GradientDescent')
        self.potential = potential
        self.loss = self.potential.loss
        self.dx = self.potential.x1
        self.train_step = self.minimize(self.loss)
        self.nfev = tf.Variable(tf.constant(0), tf.int32)
        self.rms = tf.Variable(tf.constant(1), tf.float64)
        self.ndim = tf.constant(self.potential.nspins, tf.float64)
        self.tol = tf.constant(tol, tf.float64)
        self.success = tf.constant(False, tf.bool)
        self.session = self.potential.session
        init = tf.initialize_all_variables()
        self.session.run(init)

    def _valid_dtypes(self):
        return set([tf.float32, tf.float64])

    def process_grad(self, grad, coords):
        """
        only grad can be modified because coords needs to remain a variable
        if modified it's no longer a variable
        :param grad:
        :param coords:
        :return:
        """
        zerov = tf.mul(coords, rnorm(coords))
        grad = orthog(grad, zerov)
        return grad, coords

    def compute_rms(self, grad):
        return tf.sqrt(tf.div(tf.reduce_sum(tf.square(grad)), self.ndim))

    def set_coords(self, coords):
        self.session.run(self.dx.assign(coords))

    def stop_criterion_satisfied(self, grad):
        self.rms = self.compute_rms(grad)
        return tf.less_equal(self.rms, self.tol)

    def one_iteration(self):
        grad, coords = self.compute_gradients(self.loss, var_list=[self.dx])[0]
        grad, coords = self.process_grad(grad, coords)
        self.success = self.stop_criterion_satisfied(grad)
        return self.apply_gradients([(grad, coords)], global_step=self.nfev)

    def run(self, coords, niter=100):
        self.set_coords(coords)
        for i in xrange(niter):
            self.session.run(tf.group(self.dx.assign(tf.mul(tf.mul(self.dx, rnorm(self.dx)), self.ndim)),
                                      self.one_iteration()))
            print self.rms.eval(session=self.session)
            if self.success.eval(session=self.session):
                break
        print i

    def get_results(self):
        res = Result()
        with self.session.as_default():
            g, x = self.compute_gradients(self.loss, var_list=[self.dx])[0]
            res.energy = self.loss.eval()
            res.coords = x.eval()
            res.grad = g.eval()
            res.nfev = self.nfev.eval   ()
            res.rms = self.rms.eval()
            #res.success = bool(self.thisptr.get().success())
            #res.nsteps = self.thisptr.get().get_niter()
        return res

if __name__ == "__main__":
    n=30
    p=4
    # interactions = np.ones((n,n,n))

    dtype = 'float64'
    tf.set_random_seed(100)
    norm = tf.random_normal([n for _ in xrange(p)], mean=0, stddev=1.0, dtype=dtype)
    sess = tf.Session()
    with sess.as_default():
        interactions = norm.eval()
    for comb in combinations(range(n), p):
        w = interactions[comb]
        for perm in permutations(comb):
            interactions[perm] = w

    coords = np.random.rand(n)
    coords /= np.linalg.norm(coords) / coords.size

    potTF = MeanFieldPSpinSphericalTF(interactions, n, p, dtype=dtype)
    potPL = MeanFieldPSpinSpherical(factorial(p)*interactions.flatten(), n, p)

    # gd = GradientDescent(potTF)
    # gd.run(coords)
    # print gd.get_results()

    e, grad = potPL.getEnergyGradient(coords)
    print '{0:.15f} {0:.15f}'.format(e, np.linalg.norm(grad))
    e, grad = potTF.getEnergyGradient(coords)
    print '{0:.15f} {0:.15f}'.format(e, np.linalg.norm(grad))

    import timeit
    print timeit.timeit('e, grad = potPL.getEnergyGradient(coords)', "from __main__ import potPL, coords", number=10)
    print timeit.timeit('e, grad = potTF.getEnergyGradient(coords)', "from __main__ import potTF, coords", number=10)

    def minimize(pot, coords):
       #print coords
       #print "start energy", pot.getEnergy(coords)
       results = lbfgs_cpp(coords, pot, nsteps=1e5, tol=1e-9, iprint=-1, maxstep=100)
       #results = modifiedfire_cpp(coords, pot, nsteps=1e5, tol=1e-5, iprint=-1)
       #print "quenched energy", results.energy
       results.coords /= (np.linalg.norm(results.coords)/np.sqrt(n))
       print "E: {}, nsteps: {}".format(results.energy, results.nfev)
       if results.success:
           return [results.coords, results.energy, results.nfev]

    print minimize(potPL, coords)
    print minimize(potTF, coords)
    #
    # print timeit.timeit('minimize(potPL, coords)', "from __main__ import potPL, coords, minimize", number=1)
    # print timeit.timeit('minimize(potTF, coords)', "from __main__ import potTF, coords, minimize", number=1)
