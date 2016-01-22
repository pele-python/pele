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
from pele.optimize._quench import lbfgs_cpp, modifiedfire_cpp
from pele.optimize import Result
from tensorflow.examples.tutorials.mnist import input_data

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
        grad = tf.gradients(self.loss, self.x1, colocate_gradients_with_ops=True)[0]
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

class MnistSoftmaxRegressionTF(BasePotential):
    """
    784 is the dimensionality of a flattened mnist array
    we have 784 input features and 10 outputs
    """

    def __init__(self, batchsize=3000):
        self.dtype = 'float64'
        self.mnist = input_data.read_data_sets('MNIST_data', one_hot=True)
        self.batchsize = batchsize
        self.featdim = 784
        self.outdim = 10
        self.x = tf.constant(self.mnist.train.images[:batchsize], tf.float32)
        self.y_ = tf.constant(self.mnist.train.labels[:batchsize], tf.float32)
        self.W = tf.Variable(tf.zeros([784,10])) #parameters
        self.b = tf.Variable(tf.zeros([10])) #bias
        self.session = tf.Session()
        init = tf.initialize_all_variables()
        self.session.run(init)

    def __exit__(self, exc_type, exc_value, traceback):
        self.session.close()

    @property
    def model(self):
        """
        softmax regression
        """
        return tf.matmul(self.x,self.W) + self.b

    @property
    def loss(self):
        """
        cross entropy
        """
        return tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(self.model, self.y_))

    @property
    def compute_gradient_W(self):
        return tf.gradients(self.loss, [self.W], colocate_gradients_with_ops=True)[0]

    @property
    def compute_gradient_b(self):
        return tf.gradients(self.loss, [self.b], colocate_gradients_with_ops=True)[0]

    @jit
    def pele_to_tf(self, coords):
        W = coords[:self.featdim*self.outdim].reshape((self.featdim, self.outdim))
        b = coords[-self.outdim:]
        return W, b

    @jit
    def tf_to_pele(self, W, b):
        return np.append(W.flatten(), b)

    def getEnergy(self, coords):
        W, b = self.pele_to_tf(coords)
        return self.session.run(self.loss, feed_dict={self.W : W,
                                                      self.b : b})

    def getEnergyGradient(self, coords):
        W, b = self.pele_to_tf(coords)
        e, gradW, gradb = self.session.run([self.loss, self.compute_gradient_W, self.compute_gradient_b],
                                           feed_dict={self.W : W,
                                                      self.b : b})
        return e, self.tf_to_pele(gradW, gradb)

#Build a Multilayer Convolutional Network

def weight_variable(shape):
    initial = tf.truncated_normal(shape, stddev=0.1)
    return tf.Variable(initial)

def bias_variable(shape):
    initial = tf.constant(0.1, shape=shape)
    return tf.Variable(initial)

def conv2d(x, W):
    return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

def max_pool_2x2(x):
    return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],
                          strides=[1, 2, 2, 1], padding='SAME')


class MnistMCNTF(BasePotential):
    """
    Deep convolutional network model from
    https://www.tensorflow.org/versions/master/tutorials/mnist/pros/index.html#train-the-model
    """

    def __init__(self, batchsize=50):
        self.dtype = 'float64'
        self.mnist = input_data.read_data_sets('MNIST_data', one_hot=True)
        self.batchsize = batchsize
        self.x = tf.constant(self.mnist.train.images[:batchsize], tf.float32)
        self.y_ = tf.constant(self.mnist.train.labels[:batchsize], tf.float32)
        #First Convolutional Layer
        self.W_conv1 = weight_variable([5, 5, 1, 32])
        self.b_conv1 = bias_variable([32])
        #Second Convolutional Layer
        self.W_conv2 = weight_variable([5, 5, 32, 64])
        self.b_conv2 = bias_variable([64])
        #Densely Connected Layer
        self.W_fc1 = weight_variable([7 * 7 * 64, 1024])
        self.b_fc1 = bias_variable([1024])
        #Dropout
        self.keep_prob = tf.placeholder("float")
        #Readout Layer
        self.W_fc2 = weight_variable([1024, 10])
        self.b_fc2 = bias_variable([10])

        self.session = tf.Session()
        self.train_step = tf.train.AdagradOptimizer(1e-4).minimize(self.loss)
        init = tf.initialize_all_variables()
        self.session.run(init)

    def __exit__(self, exc_type, exc_value, traceback):
        self.session.close()

    @property
    def model(self):
        """
        Multilayer Convolutional Network
        """
        x_image = tf.reshape(self.x, [-1,28,28,1])
        h_conv1 = tf.nn.relu(conv2d(x_image, self.W_conv1) + self.b_conv1)
        h_pool1 = max_pool_2x2(h_conv1)
        h_conv2 = tf.nn.relu(conv2d(h_pool1, self.W_conv2) + self.b_conv2)
        h_pool2 = max_pool_2x2(h_conv2)
        h_pool2_flat = tf.reshape(h_pool2, [-1, 7*7*64])
        h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, self.W_fc1) + self.b_fc1)
        h_fc1_drop = tf.nn.dropout(h_fc1, self.keep_prob)
        return tf.matmul(h_fc1_drop, self.W_fc2) + self.b_fc2

    @property
    def loss(self):
        """
        cross entropy
        """
        return tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(self.model, self.y_))

    def compute_gradient(self, var):
        return tf.gradients(self.loss, [var], colocate_gradients_with_ops=True)[0]

    @jit
    def pele_to_tf(self, coords):
        nconv1 = 5*5*1*32
        W_conv1 = coords[:nconv1].reshape([5, 5, 1, 32])
        b_conv1 = coords[nconv1:nconv1+32]
        #Second Convolutional Layer
        nconv2 = 5*5*1*32 + 32 + 5*5*32*64
        W_conv2 = coords[nconv1+32:nconv2].reshape([5, 5, 32, 64])
        b_conv2 = coords[nconv2:nconv2+64]
        #Densely Connected Layer
        nconv3 = nconv2 + 64 + 7 * 7 * 64 * 1024
        W_fc1 = coords[nconv2+64:nconv3].reshape([7 * 7 * 64, 1024])
        b_fc1 = coords[nconv3:nconv3+1024]
        #Readout Layer
        nconv4 = nconv3+1024 + 1024 * 10
        W_fc2 = coords[nconv3+1024:nconv4].reshape([1024, 10])
        b_fc2 = coords[nconv4:]
        vars = [W_conv1, b_conv1, W_conv2, b_conv2, W_fc1, b_fc1, W_fc2, b_fc2]
        return vars

    @jit
    def tf_to_pele(self, vars):
        coords = np.empty(0)
        for var in vars:
            coords = np.append(coords.flatten(), var.flatten())
        return coords

    def getEnergy(self, coords):
        vars = self.pele_to_tf(coords)
        #gradop = [self.compute_gradient(var)[0] for var in vars]
        return self.session.run(self.loss,
                               feed_dict={self.W_conv1 : vars[0],
                                          self.b_conv1 : vars[1],
                                          self.W_conv2 : vars[2],
                                          self.b_conv2 : vars[3],
                                          self.W_fc1 : vars[4],
                                          self.b_fc1 : vars[5],
                                          self.W_fc2 : vars[6],
                                          self.b_fc2 : vars[7],
                                          self.keep_prob: 1})

    def getEnergyGradient(self, coords):
        vars = self.pele_to_tf(coords)
        out = self.session.run([self.loss,
                                self.compute_gradient(self.W_conv1),
                                self.compute_gradient(self.b_conv1),
                                self.compute_gradient(self.W_conv2),
                                self.compute_gradient(self.b_conv2),
                                self.compute_gradient(self.W_fc1),
                                self.compute_gradient(self.b_fc1),
                                self.compute_gradient(self.W_fc2),
                                self.compute_gradient(self.b_fc2)
                                ],
                               feed_dict={self.W_conv1 : vars[0],
                                          self.b_conv1 : vars[1],
                                          self.W_conv2 : vars[2],
                                          self.b_conv2 : vars[3],
                                          self.W_fc1 : vars[4],
                                          self.b_fc1 : vars[5],
                                          self.W_fc2 : vars[6],
                                          self.b_fc2 : vars[7],
                                          self.keep_prob: 1})
        return out[0], self.tf_to_pele(out[1:])

    def minimize(self):
        for i in range(100):
            self.session.run(self.train_step, feed_dict={self.keep_prob : 0.5})
        print self.session.run(self.loss, feed_dict={self.keep_prob : 0.5})

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
    import timeit
    np.random.seed(43)

    def minimize(pot, coords):
        print "start energy", pot.getEnergy(coords)
        results = lbfgs_cpp(coords, pot, M=4, nsteps=1e5, tol=1e-5, iprint=1, verbosity=5, maxstep=1)
        #results = modifiedfire_cpp(coords, pot, nsteps=1e5, tol=1e-5, iprint=1, verbosity=5, maxstep=1)
        print "quenched energy", results.energy
        print "E: {}, nsteps: {}".format(results.energy, results.nfev)
        if results.success:
           return [results.coords, results.energy, results.nfev]

    if False:
        pot = MnistMCNTF()
        pot.minimize()

    if False:
        pot = MnistMCNTF()
        W_conv1 = np.abs(np.random.normal(scale=0.1, size=[5, 5, 1, 32]))
        b_conv1 = np.ones([32]) * 0.1
        #Second Convolutional Layer
        W_conv2 = np.abs(np.random.normal(scale=0.1, size=[5, 5, 32, 64]))
        b_conv2 = np.ones([64]) * 0.1
        #Densely Connected Layer
        W_fc1 = np.abs(np.random.normal(scale=0.1, size=[7 * 7 * 64, 1024]))
        b_fc1 = np.ones([1024]) * 0.1
        #Readout Layer
        W_fc2 = np.abs(np.random.normal(scale=0.1, size=[1024, 10]))
        b_fc2 = np.ones([10]) * 0.1
        vars = [W_conv1, b_conv1, W_conv2, b_conv2, W_fc1, b_fc1, W_fc2, b_fc2]
        coords = np.empty(0)
        for var in vars:
            coords = np.append(coords.flatten(), var.flatten())
        print coords.shape
        print pot.getEnergy(coords)
        print pot.getEnergyGradient(coords)
        minimize(pot, coords)

    if True:
        pot = MnistSoftmaxRegressionTF()
        w = np.abs(np.random.normal(size=[784,10], scale=0.1).flatten())
        b = np.ones(10)*0.1
        coords = np.append(w, b)
        #print pot.getEnergy(coords)
        #print pot.getEnergyGradient(coords)
        minimize(pot, coords)

    if False:
        dtype = 'float64'
        n=100
        p=3
        # interactions = np.ones((n,n,n))
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

        # print timeit.timeit('e, grad = potPL.getEnergyGradient(coords)', "from __main__ import potPL, coords", number=10)
        # print timeit.timeit('e, grad = potTF.getEnergyGradient(coords)', "from __main__ import potTF, coords", number=10)

        #print minimize(potPL, coords)
        print minimize(potTF, coords)
        #
        # print timeit.timeit('minimize(potPL, coords)', "from __main__ import potPL, coords, minimize", number=1)
        #print timeit.timeit('minimize(potTF, coords)', "from __main__ import potTF, coords, minimize", number=1)
