from past.utils import old_div
import numpy

# LEPS 2d potential
class leps:
    def getEnergy(self, r):
        """
        potential energy as a function of position
        for the LEPS potential on a line
        python version
        """
        x = r[0]
        y = r[1]
        a = 0.05
        b = 0.3
        c = 0.05
        alpha = 1.942
        r0 = 0.742
        dAB = 4.746
        dBC = 4.746
        dAC = 3.445

        def Q(d, r):
            return old_div(
                d
                * (
                    old_div(3 * numpy.exp(-2 * alpha * (r - r0)), 2)
                    - numpy.exp(-alpha * (r - r0))
                ),
                2,
            )

        def J(d, r):
            return old_div(
                d
                * (
                    numpy.exp(-2 * alpha * (r - r0))
                    - 6 * numpy.exp(-alpha * (r - r0))
                ),
                4,
            )

        rAB = x
        rBC = y
        rAC = rAB + rBC

        JABred = old_div(J(dAB, rAB), (1 + a))
        JBCred = old_div(J(dBC, rBC), (1 + b))
        JACred = old_div(J(dAC, rAC), (1 + c))

        return (
            old_div(Q(dAB, rAB), (1 + a))
            + old_div(Q(dBC, rBC), (1 + b))
            + old_div(Q(dAC, rAC), (1 + c))
            - numpy.sqrt(
                JABred * JABred
                + JBCred * JBCred
                + JACred * JACred
                - JABred * JBCred
                - JBCred * JACred
                - JABred * JACred
            )
        )

    def getEnergyGradient(self, r):
        """
        force as a function of position
        for the LEPS potential on a line
        python version
        """
        x = r[0]
        y = r[1]
        a = 0.05
        b = 0.3
        c = 0.05
        alpha = 1.942
        r0 = 0.742
        dAB = 4.746
        dBC = 4.746
        dAC = 3.445

        def Q(d, r):
            return old_div(
                d
                * (
                    old_div(3 * numpy.exp(-2 * alpha * (r - r0)), 2)
                    - numpy.exp(-alpha * (r - r0))
                ),
                2,
            )

        def J(d, r):
            return old_div(
                d
                * (
                    numpy.exp(-2 * alpha * (r - r0))
                    - 6 * numpy.exp(-alpha * (r - r0))
                ),
                4,
            )

        def dQ(d, r):
            return old_div(
                alpha
                * d
                * (
                    -3 * numpy.exp(-2 * alpha * (r - r0))
                    + numpy.exp(-alpha * (r - r0))
                ),
                2,
            )

        def dJ(d, r):
            return old_div(
                alpha
                * d
                * (
                    -2 * numpy.exp(-2 * alpha * (r - r0))
                    + 6 * numpy.exp(-alpha * (r - r0))
                ),
                4,
            )

        rAB = x
        rBC = y
        rAC = rAB + rBC

        JABred = old_div(J(dAB, rAB), (1 + a))
        JBCred = old_div(J(dBC, rBC), (1 + b))
        JACred = old_div(J(dAC, rAC), (1 + c))

        dJABred = old_div(dJ(dAB, rAB), (1 + a))
        dJBCred = old_div(dJ(dBC, rBC), (1 + b))
        dJACred = old_div(dJ(dAC, rAC), (1 + c))

        Fx = (
            old_div(dQ(dAB, rAB), (1 + a))
            + old_div(dQ(dAC, rAC), (1 + c))
            - old_div(
                (
                    2 * JABred * dJABred
                    + 2 * JACred * dJACred
                    - dJABred * JBCred
                    - JBCred * dJACred
                    - dJABred * JACred
                    - JABred * dJACred
                ),
                (
                    2
                    * numpy.sqrt(
                        JABred * JABred
                        + JBCred * JBCred
                        + JACred * JACred
                        - JABred * JBCred
                        - JBCred * JACred
                        - JABred * JACred
                    )
                ),
            )
        )

        Fy = (
            old_div(dQ(dBC, rBC), (1 + b))
            + old_div(dQ(dAC, rAC), (1 + c))
            - old_div(
                (
                    2 * JBCred * dJBCred
                    + 2 * JACred * dJACred
                    - JABred * dJBCred
                    - dJBCred * JACred
                    - JBCred * dJACred
                    - JABred * dJACred
                ),
                (
                    2
                    * numpy.sqrt(
                        JABred * JABred
                        + JBCred * JBCred
                        + JACred * JACred
                        - JABred * JBCred
                        - JBCred * JACred
                        - JABred * JACred
                    )
                ),
            )
        )

        return self.getEnergy(r), numpy.array([Fx, Fy])
