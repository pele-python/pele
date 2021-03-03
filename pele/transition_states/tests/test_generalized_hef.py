from __future__ import absolute_import
import unittest

from pele.transition_states._generalized_dimer import GeneralizedDimer
from . import test_generalized_dimer


class TestGeneralizedHEF(test_generalized_dimer.TestGeneralizedDimer):
    def make_dimer(self, x):
        return GeneralizedDimer(x, self.pot,
                                leig_kwargs=dict(orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors(),
                                                 iprint=10,
                                                 ) ,
                                translator_kwargs=dict(iprint=1
                                                       ),
                                dimer=False,
#                                maxiter=100,
                                translational_steps=3,
                                )
    

if __name__ == "__main__":
    unittest.main()

