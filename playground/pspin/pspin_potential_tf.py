import numpy as np
import tensorflow as tf
from pele.potentials import BasePotential

class DoubleGDOptimizer(tf.train.GradientDescentOptimizer):
    """
    this wrapper adds support for float64 for tf.train.GradientDescentOptimizer
    """
    def _valid_dtypes(self):
        return set([tf.float32, tf.float64])

class MeanFieldPSpinSphericalTF(BasePotential):
    def __init__(self, nspins, interaction, p=3):
        assert p==3
        assert interaction.shape == (nspins, nspins, nspins)
        self.nspins = nspins
        self.interaction = tf.constant(interaction, tf.float64)
        # self.x_ = tf.placeholder(dtype='float64', shape=[nspins])
        self.xi = tf.Variable(tf.zeros([self.nspins], dtype='float64'))
        self.xj = tf.Variable(tf.zeros([self.nspins], dtype='float64'))
        self.xk = tf.Variable(tf.zeros([self.nspins], dtype='float64'))
        self.session = tf.Session()
        init = tf.initialize_all_variables()
        self.session.run(init)

    def __exit__(self, exc_type, exc_value, traceback):
        self.session.close()

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
        return tf.reduce_sum(tf.mul(self.interaction, self.lossTensor))

    def getEnergy(self, coords):
        self.session.run(tf.group(self.xi.assign(coords),
                                  self.xj.assign(coords),
                                  self.xk.assign(coords)))
        e = self.session.run(self.loss)
        return -e

    def getEnergyGradient(self, coords):
        self.session.run(tf.group(self.xi.assign(coords),
                                  self.xj.assign(coords),
                                  self.xk.assign(coords)))
        e = self.session.run(self.loss)
        print tf.gradients(self.loss, self.xk)[0].eval(session=self.session)
        print self.xk.eval(session=self.session)
        return -e


if __name__ == "__main__":
    n = 5
    pot = MeanFieldPSpinSphericalTF(n, np.ones([n,n,n], dtype='float64'))
    coords = np.ones(n, dtype='float64')
    print pot.getEnergy(coords)
    print pot.getEnergyGradient(coords)