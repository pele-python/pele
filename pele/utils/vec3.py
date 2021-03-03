import numpy as np

__all__ = ["invert3x3"]


def invert3x3(A):
    """return the inverse of a 3x3 matrix
    """
    Ainv = np.zeros([3, 3])
    Ainv[0, 0] = A[1, 1] * A[2, 2] - A[2, 1] * A[1, 2]
    Ainv[0, 1] = A[2, 1] * A[0, 2] - A[0, 1] * A[2, 2]
    Ainv[0, 2] = A[0, 1] * A[1, 2] - A[0, 2] * A[1, 1]

    Adet = Ainv[0, 0] * A[0, 0] + Ainv[0, 1] * A[1, 0] + Ainv[0, 2] * A[2, 0]

    Ainv[0, 0] = Ainv[0, 0] / Adet
    Ainv[0, 1] = Ainv[0, 1] / Adet
    Ainv[0, 2] = Ainv[0, 2] / Adet

    Ainv[1, 0] = (A[1, 2] * A[2, 0] - A[1, 0] * A[2, 2]) / Adet
    Ainv[1, 1] = (A[0, 0] * A[2, 2] - A[2, 0] * A[0, 2]) / Adet
    Ainv[1, 2] = (A[1, 0] * A[0, 2] - A[0, 0] * A[1, 2]) / Adet

    Ainv[2, 0] = (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]) / Adet
    Ainv[2, 1] = (A[2, 0] * A[0, 1] - A[0, 0] * A[2, 1]) / Adet
    Ainv[2, 2] = (A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]) / Adet
    return Ainv

