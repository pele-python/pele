import numpy as np
import tensorflow as tf
from pele.potentials import BasePotential

class MeanFieldPSpinSphericalTF(BasePotential):
    def __init__(self, nspins, interaction, p=3):
        assert p==3
        assert interaction.shape == (nspins, nspins, nspins)
        self.nspins = nspins
        self.interaction = tf.constant(interaction)
        self.x_ = tf.placeholder(dtype='float32', shape=[nspins])
        #self.x_ = tf.Variable(tf.zeros([self.nspins], dtype='float32'))
        init = tf.initialize_all_variables()
        self.session = tf.Session()
        self.session.run(init)

    def __exit__(self, exc_type, exc_value, traceback):
        self.session.close()

    @property
    def lossTensorPartial(self):
        """this gives a 2x2 matrix generator"""
        return tf.mul(tf.reshape(self.x_, [self.nspins]), tf.reshape(self.x_, [self.nspins,1]))

    @property
    def lossTensor(self):
        """this gives the full loss tensor generator"""
        return tf.mul(self.lossTensorPartial, tf.reshape(self.x_, [self.nspins,1,1]))

    @property
    def loss(self):
        """this mutliplies times the interaction and reduces sum reduces the tensor"""
        return tf.reduce_sum(tf.mul(self.interaction, self.lossTensor))

    def getEnergy(self, coords):
        e = self.session.run(self.loss, feed_dict={self.x_:coords})
        return -e

    # @property
    # def loss_grad(self):
    #     return tf.gradients(self.lossTensor, self.x_)

    # def getEnergyGradient(self, coords):
    #   e = self.session.run(self.loss, feed_dict={self.x_:coords})
    #   g = self.session.run(self.loss_grad, feed_dict={self.x_:coords})
    #   grads_and_vars = self.opt.compute_gradients(self.loss)
    #   self.opt = tf.train.GradientDescentOptimizer(learning_rate=0.1)
    #   g, v = self.session.run(grads_and_vars)
    #   for g, v in grads_and_vars:
    #       print v
    #       print g

if __name__ == "__main__":
    pot = MeanFieldPSpinSphericalTF(5, np.ones([5,5,5], dtype='float32'))
    e = pot.getEnergy(np.ones(5, dtype='float32'))
    print e
    assert e == -125.