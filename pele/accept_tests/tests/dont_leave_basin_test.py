import unittest
import random
import pele.exceptions as exc
import pele.accept_tests.dont_leave_basin as dlb


class TestDontLeaveBasinTest(unittest.TestCase):
    def setUp(self):
        # create some DontLeaveBasin objects with different energy
        # criteria
        self.dlb_e_plus2 = dlb.DontLeaveBasin(Ecriterion=100.0)
        self.dlb_e_0 = dlb.DontLeaveBasin(Ecriterion=1.0)
        self.dlb_e_4 = dlb.DontLeaveBasin(Ecriterion=1.0e-4)
        self.dlb_e_6 = dlb.DontLeaveBasin(Ecriterion=1.0e-6)
        self.dlb_default = dlb.DontLeaveBasin()

    def test_default(self):
        # test that the default DontLeaveBasin object has the same Ecriterion
        # as dlb_e_4
        self.assertEquals(self.dlb_e_4.Ecriterion,
                          self.dlb_default.Ecriterion)

    def test_acceptPosPos(self):
        # test that cases when both energies are positive and the difference is
        # under the threshold are accepted
        Enew = random.uniform(100, 1000.0)
        Eold = Enew - 100.0 * random.random()
        self.assertTrue(self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_plus2.acceptReject(Enew, Eold))
        Eold = Enew - 1.0 * random.random()
        self.assertTrue(self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_0.acceptReject(Enew, Eold))
        Eold = Enew - 1.0e-4 * random.random()
        self.assertTrue(self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_4.acceptReject(Enew, Eold))
        Eold = Enew - 1.0e-6 * random.random()
        self.assertTrue(self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_6.acceptReject(Enew, Eold))

    def test_acceptNegNeg(self):
        # test that cases when both energies are negative and the difference is
        # under the threshold are accepted
        Enew = -random.uniform(100, 1000.0)
        Eold = Enew + 100.0 * random.random()
        self.assertTrue(self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_plus2.acceptReject(Enew, Eold))
        Eold = Enew + 1.0 * random.random()
        self.assertTrue(self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_0.acceptReject(Enew, Eold))
        Eold = Enew + 1.0e-4 * random.random()
        self.assertTrue(self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_4.acceptReject(Enew, Eold))
        Eold = Enew + 1.0e-6 * random.random()
        self.assertTrue(self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_6.acceptReject(Enew, Eold))

    def test_acceptPosNeg(self):
        # test that cases when one energy is positive and the other is negative
        # and the difference is under the threshold are accepted
        Enew = random.uniform(0, 50)
        Eold = random.uniform(-50, 0)
        self.assertTrue(self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_plus2.acceptReject(Enew, Eold))
        Enew = random.uniform(0, 0.5)
        Eold = random.uniform(-0.5, 0)
        self.assertTrue(self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_0.acceptReject(Enew, Eold))
        Enew = random.uniform(0, 0.5e-4)
        Eold = random.uniform(-0.5e-4, 0)
        self.assertTrue(self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_4.acceptReject(Enew, Eold))
        Enew = random.uniform(0, 0.5e-6)
        Eold = random.uniform(-0.5e-6, 0)
        self.assertTrue(self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_6.acceptReject(Enew, Eold))

    def test_acceptZeroPos(self):
        # test that cases when one energy is zero and the other is positive 
        # and under the threshold are accepted
        Eold = 0
        Enew = 100.0 * random.random()
        self.assertTrue(self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_plus2.acceptReject(Enew, Eold))
        Enew = 1.0 * random.random()
        self.assertTrue(self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_0.acceptReject(Enew, Eold))
        Enew = 1.0e-4 * random.random()
        self.assertTrue(self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_4.acceptReject(Enew, Eold))
        Enew = 1.0e-6 * random.random()
        self.assertTrue(self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_6.acceptReject(Enew, Eold))

    def test_acceptZeroNeg(self):
        # test that cases when one energy is zero and the other is negative 
        # and under the threshold are accepted
        Eold = 0
        Enew = -100.0 * random.random()
        self.assertTrue(self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_plus2.acceptReject(Enew, Eold))
        Enew = -1.0 * random.random()
        self.assertTrue(self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_0.acceptReject(Enew, Eold))
        Enew = -1.0e-4 * random.random()
        self.assertTrue(self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_4.acceptReject(Enew, Eold))
        Enew = -1.0e-6 * random.random()
        self.assertTrue(self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(self.dlb_e_6.acceptReject(Enew, Eold))

    def test_rejectPosPos(self):
        # test that cases when both energies are positive and the difference is
        # under the threshold are accepted
        Enew = random.uniform(100, 1000.0)
        Eold = Enew - 100.0 / random.random()
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Enew, Eold))
        Eold = Enew - 1.0 / random.random()
        self.assertTrue(not self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_0.acceptReject(Enew, Eold))
        Eold = Enew - 1.0e-4 / random.random()
        self.assertTrue(not self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_4.acceptReject(Enew, Eold))
        Eold = Enew - 1.0e-6 / random.random()
        self.assertTrue(not self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_6.acceptReject(Enew, Eold))

    def test_rejectNegNeg(self):
        # test that cases when both energies are negative and the difference is
        # under the threshold are accepted
        Enew = -random.uniform(100, 1000.0)
        Eold = Enew + 100.0 / random.random()
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Enew, Eold))
        Eold = Enew + 1.0 / random.random()
        self.assertTrue(not self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_0.acceptReject(Enew, Eold))
        Eold = Enew + 1.0e-4 / random.random()
        self.assertTrue(not self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_4.acceptReject(Enew, Eold))
        Eold = Enew + 1.0e-6 / random.random()
        self.assertTrue(not self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_6.acceptReject(Enew, Eold))

    def test_rejectPosNeg(self):
        # test that cases when one energy is positive and the other is negative
        # and the difference is under the threshold are accepted
        Enew = random.uniform(50, 100)
        Eold = random.uniform(-100, -50)
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Enew, Eold))
        Enew = random.uniform(0.5, 1.0)
        Eold = random.uniform(-1.0, -0.5)
        self.assertTrue(not self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_0.acceptReject(Enew, Eold))
        Enew = random.uniform(0.5e-4, 1.0e-4)
        Eold = random.uniform(-1.0e-4, -0.5e-4)
        self.assertTrue(not self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_4.acceptReject(Enew, Eold))
        Enew = random.uniform(0.5e-6, 1.0e-6)
        Eold = random.uniform(-1.0e-6, -0.5e-6)
        self.assertTrue(not self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_6.acceptReject(Enew, Eold))

    def test_rejectZeroPos(self):
        # test that cases when one energy is zero and the other is positive 
        # and over the threshold are rejected
        Eold = 0
        Enew = 100.0 / random.random()
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Enew, Eold))
        Enew = 1.0 / random.random()
        self.assertTrue(not self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_0.acceptReject(Enew, Eold))
        Enew = 1.0e-4 / random.random()
        self.assertTrue(not self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_4.acceptReject(Enew, Eold))
        Enew = 1.0e-6 / random.random()
        self.assertTrue(not self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_6.acceptReject(Enew, Eold))

    def test_rejectZeroNeg(self):
        # test that cases when one energy is zero and the other is negative 
        # and over the threshold are rejected
        Eold = 0
        Enew = -100.0 / random.random()
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_plus2.acceptReject(Enew, Eold))
        Enew = -1.0 / random.random()
        self.assertTrue(not self.dlb_e_0.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_0.acceptReject(Enew, Eold))
        Enew = -1.0e-4 / random.random()
        self.assertTrue(not self.dlb_e_4.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_4.acceptReject(Enew, Eold))
        Enew = -1.0e-6 / random.random()
        self.assertTrue(not self.dlb_e_6.acceptReject(Eold, Enew))
        self.assertTrue(not self.dlb_e_6.acceptReject(Enew, Eold))

    def test_rejectZeroCriterion(self):
        # test that with an energy criterion of zero, all changes in energy
        # are rejected
        # Create a DontLeaveBasin object with Ecriterion = 0
        dlb_zero = dlb.DontLeaveBasin(Ecriterion=0.0)
        (Eold, Enew) = (random.random(), random.random())
        self.assertTrue(not dlb_zero.acceptReject(Eold, Enew))

    def test_negativeCriterionException(self):
        # test that trying to use a negative energy criterion raises the 
        # SignError exception 
        self.assertRaises(exc.SignError, dlb.DontLeaveBasin, -100.0)
        self.assertRaises(exc.SignError, dlb.DontLeaveBasin, -1.0)
        self.assertRaises(exc.SignError, dlb.DontLeaveBasin, -1.0e-4)
        self.assertRaises(exc.SignError, dlb.DontLeaveBasin, -1.0e-6)


if __name__ == '__main__':
    unittest.main()
