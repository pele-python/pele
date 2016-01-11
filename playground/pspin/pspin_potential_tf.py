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

class DoubleGDOptimizer(tf.train.GradientDescentOptimizer):
    """
    this wrapper adds support for float64 for tf.train.GradientDescentOptimizer
    """
    def _valid_dtypes(self):
        return set([tf.float32, tf.float64])

def MeanFieldPSpinSphericalTF(interactions, nspins, p):
    if p == 3:
        return MeanField3SpinSphericalTF(interactions, nspins)
    elif p == 4:
        return MeanField4SpinSphericalTF(interactions, nspins)
    elif p == 5:
        return MeanField5SpinSphericalTF(interactions, nspins)
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
        self.prf =  np.power(np.sqrt(self.nspins), self.p-1) if self.p > 2 else 1.
        self.sqrtN = np.sqrt(self.nspins)
        interactions = self._adaptInteractions(interactions)
        # print len(np.flatnonzero(interactions[0]))
        self.interactions = tf.constant(interactions, dtype=self.dtype)
        self.x1 = tf.Variable(tf.zeros([self.nspins], dtype=self.dtype))
        self.x2 = tf.Variable(tf.zeros([self.nspins], dtype=self.dtype))
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
        return tf.reduce_sum(tf.mul(self.interactions, self.lossTensor))

    # @jit
    def _orthog_to_zero(self, v, zerov):
        return orthogonalize(v, [zerov / np.linalg.norm(zerov)])

    # @jit
    def _normalizeSpins(self, coords):
        coords /= (np.linalg.norm(coords)/self.sqrtN)

    def _setSpinsTf(self, coords):
        self.session.run(tf.group(self.x1.assign(coords),
                                  self.x2.assign(coords)))

    def _setSpins(self, coords):
        self._normalizeSpins(coords)
        self._setSpinsTf(coords)

    def getEnergy(self, coords):
        self._setSpins(coords)
        e = -self.session.run(self.loss)
        return e/self.prf

    def getEnergyGradient(self, coords):
        self._setSpins(coords)
        e = -self.session.run(self.loss)
        grad = -tf.gradients(self.loss, self.x1)[0].eval(session=self.session)
        grad = self._orthog_to_zero(grad/self.prf, coords)
        return e/self.prf, grad

class MeanField3SpinSphericalTF(BaseMeanFieldPSpinSphericalTF):
    def __init__(self, interactions, nspins):
        super(MeanField3SpinSphericalTF, self).__init__(interactions, nspins, p=3)

    @property
    def lossTensor(self):
        """this gives the full loss tensor generator"""
        return tf.mul(self.lossTensorPartial, tf.reshape(self.x2, [self.nspins,1,1]))

class MeanField4SpinSphericalTF(BaseMeanFieldPSpinSphericalTF):
    def __init__(self, interactions, nspins):
        super(MeanField4SpinSphericalTF, self).__init__(interactions, nspins, p=4)

    @property
    def lossTensor(self):
        """this gives the full loss tensor generator"""
        return tf.mul(tf.mul(self.lossTensorPartial, tf.reshape(self.x2, [self.nspins,1,1])),
                      tf.reshape(self.x2, [self.nspins,1,1,1]))

class MeanField5SpinSphericalTF(BaseMeanFieldPSpinSphericalTF):
    def __init__(self, interactions, nspins):
        super(MeanField5SpinSphericalTF, self).__init__(interactions, nspins, p=5)

    @property
    def lossTensor(self):
        """this gives the full loss tensor generator"""
        return tf.mul(tf.mul(tf.mul(self.lossTensorPartial, tf.reshape(self.x2, [self.nspins,1,1])),
                             tf.reshape(self.x2, [self.nspins,1,1,1])),
                      tf.reshape(self.x2, [self.nspins,1,1,1,1]))

if __name__ == "__main__":
    n=100
    p=3
    # interactions = np.ones((n,n,n))
    np.random.seed(100)
    interactions = np.zeros([n for i in xrange(p)])
    for comb in combinations(range(n), p):
        w = np.random.normal(0, np.sqrt(factorial(p)))
        for perm in permutations(comb):
            interactions[perm] = w

    potTF = MeanFieldPSpinSphericalTF(interactions/factorial(p), n, p)
    potPL = MeanFieldPSpinSpherical(interactions.flatten(), n, p)
    coords = np.ones(n, dtype='float64')

    e, grad = potPL.getEnergyGradient(coords)
    print e, np.linalg.norm(grad)
    e, grad = potTF.getEnergyGradient(coords)
    print e, np.linalg.norm(grad)

    import timeit
    print timeit.timeit('e, grad = potPL.getEnergyGradient(coords)', "from __main__ import potPL, coords", number=20)
    print timeit.timeit('e, grad = potTF.getEnergyGradient(coords)', "from __main__ import potTF, coords", number=20)

    def minimize(pot, coords):
        #print coords
        #print "start energy", pot.getEnergy(coords)
        results = lbfgs_cpp(coords, pot, nsteps=1e5, tol=1e-9, iprint=-1, maxstep=10)
        #results = modifiedfire_cpp(coords, pot, nsteps=1e5, tol=1e-5, iprint=-1)
        #print "quenched energy", results.energy
        results.coords /= (np.linalg.norm(results.coords)/np.sqrt(n))
        if results.success:
            return [results.coords, results.energy, results.nfev]

    # print minimize(potPL, coords)
    print minimize(potTF, coords)

    # print timeit.timeit('minimize(potPL, coords)', "from __main__ import potPL, coords, minimize", number=1)
    # print timeit.timeit('minimize(potTF, coords)', "from __main__ import potTF, coords, minimize", number=5)