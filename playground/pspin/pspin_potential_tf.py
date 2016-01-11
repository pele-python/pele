import numpy as np
import tensorflow as tf
from pele.potentials import BasePotential
from itertools import permutations, combinations
from pele.transition_states._zeroev import orthogonalize
from scipy.special import factorial
from pele.potentials import MeanFieldPSpinSpherical
from numba import jit

class DoubleGDOptimizer(tf.train.GradientDescentOptimizer):
    """
    this wrapper adds support for float64 for tf.train.GradientDescentOptimizer
    """
    def _valid_dtypes(self):
        return set([tf.float32, tf.float64])

class MeanFieldPSpinSphericalTF(BasePotential):
    def __init__(self, interactions, nspins, p=3):
        assert p==3
        assert interactions.shape == (nspins, nspins, nspins)
        self.p = p
        self.nspins = nspins
        self.prf = factorial(self.p) * np.power(np.sqrt(self.nspins), self.p-1) if self.p > 2 else 1.
        interactions = self._adaptInteractions(interactions)
        self.interactions = tf.constant(interactions, tf.float64)
        self.xi = tf.Variable(tf.zeros([self.nspins], dtype='float64'))
        self.xj = tf.Variable(tf.zeros([self.nspins], dtype='float64'))
        self.xk = tf.Variable(tf.zeros([self.nspins], dtype='float64'))
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
        return tf.mul(tf.reshape(self.xi, [self.nspins]), tf.reshape(self.xj, [self.nspins,1]))

    @property
    def lossTensor(self):
        """this gives the full loss tensor generator"""
        return tf.mul(self.lossTensorPartial, tf.reshape(self.xk, [self.nspins,1,1]))

    @property
    def loss(self):
        """this mutliplies times the interaction and reduces sum reduces the tensor"""
        return tf.reduce_sum(tf.mul(self.interactions, self.lossTensor))

    def _setSpins(self, coords):
        self.session.run(tf.group(self.xi.assign(coords),
                                  self.xj.assign(coords),
                                  self.xk.assign(coords)))
    @jit
    def _orthog_to_zero(self, v, zerov):
        return orthogonalize(v, [zerov])

    @jit
    def getEnergy(self, coords):
        # coords /= np.linalg.norm(coords)
        self._setSpins(coords)
        e = -self.session.run(self.loss)
        return e/self.prf

    @jit
    def getEnergyGradient(self, coords):
        # coords /= np.linalg.norm(coords)
        self._setSpins(coords)
        e = -self.session.run(self.loss)
        grad = -tf.gradients(self.loss, self.xk)[0].eval(session=self.session)
        grad = self._orthog_to_zero(grad/self.prf, coords)
        return e/self.prf, grad

if __name__ == "__main__":
    n = 150
    # interactions = np.random.random((n,n,n))
    interactions = np.ones((n,n,n))

    potTF = MeanFieldPSpinSphericalTF(interactions, n)
    potPL = MeanFieldPSpinSpherical(interactions.flatten(), n, 3)
    coords = np.ones(n, dtype='float64')
    # coords /= np.linalg.norm(coords)
    import timeit
    print timeit.timeit('e, grad = potPL.getEnergyGradient(coords)', "from __main__ import potPL, coords", number=100)
    print timeit.timeit('e, grad = potTF.getEnergyGradient(coords)', "from __main__ import potTF, coords", number=100)