import unittest
import numpy as np

from pele.potentials._fin_sin import FinSin

import _base_test


_x = np.array([-1.33239375,  0.67919166, -0.87418082,  1.47156149, -1.18321309,
        0.97785102, -1.99014857,  0.39960338, -0.23665445, -1.09375111,
        0.04556779, -1.77876227, -1.13505724, -1.48824961,  1.11516898,
        0.83611026,  1.58232546,  1.00622132, -0.82877301, -0.62487541,
       -0.14168195,  1.28506775,  0.78692738, -1.17924569,  1.78384017,
        1.26595273,  1.79889982,  1.58496521, -0.5815269 , -1.66498399])


class TestFinSin(_base_test._TestConfiguration):
    def setUp(self):
        self.pot = FinSin()
        self.x0 = _x.copy()
        self.e0 = -45.02392166668808

if __name__ == "__main__":
    unittest.main()
    
