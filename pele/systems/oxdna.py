import numpy as np
from pele import takestep
from math import pi
from pele.utils.rbtools import CoordsAdapter
from pele.utils import rotations

# choose bond to displace with linear distribution from the sides of the chain.
# The end bonds are always accepted, the middle one is accepted
# with probability P_mid. In between it's linearly interpolated.
# This leads to that middle bonds are displaced less often.
# N is the number of bonds 
def choose_bond(N, P_mid=0.):
    mid = 0.5 * float(N) - 0.5

    while True:
        i = np.random.randint(N)
        dist = float(min(i, N - i - 1))
        if (1. - P_mid) * dist / mid < np.random.random():
            return i


# This is the takestep routine for OXDNA. It is a standard rigid body takestep
# routine, but I put it here to be able to start modifying it
class OXDNATakestep(takestep.TakestepInterface):
    def __init__(self, displace=1.0, rotate=0.5 * pi, rotate_around_backbone=False):
        self.displace = displace
        self.rotate = rotate
        self.rotate_around_backbone = rotate_around_backbone

    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size / 6, coords=coords)

        # random displacement for positions
        ca.posRigid[:] += 2. * self.displace * (np.random.random(ca.posRigid.shape) - 0.5)

        # determine backbone beads
        for com, p in zip(ca.posRigid, ca.rotRigid):
            p_rnd = rotations.small_random_aa(self.rotate)
            # adjust center of mass
            if self.rotate_around_backbone:
                a1 = np.dot(rotations.aa2mx(p), np.array([1., 0., 0.]))
                x1 = com - 0.4 * a1
                mx = rotations.aa2mx(p_rnd)
                com[:] = np.dot(mx, com - x1) + x1
                # random rotation for angle-axis vectors
            p[:] = rotations.rotate_aa(p, p_rnd)


    # this is necessary for adaptive step taking
    def scale(self, factor):
        self.rotate *= factor
        self.displace *= factor

    @property
    def stepsize(self):
        return [self.rotate, self.displace]


class OXDNAScrewStep(takestep.TakestepInterface):
    def __init__(self, rotate_backbone=0.5 * pi, rotate_base=0., ntorsionmoves=1, P_mid=1.):
        self.rotate_backbone = rotate_backbone
        self.rotate_base = rotate_base
        self.ntorsionmoves = ntorsionmoves
        self.P_mid = P_mid

    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size / 6, coords=coords)

        for i in range(ca.nrigid):
            a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
            a *= 2. * (np.random.random() - 0.5) * self.rotate_base
            ca.rotRigid[i] = rotations.rotate_aa(ca.rotRigid[i], a)

            # random rotation for angle-axis vectors
        if self.rotate_backbone != 0.:
            for j in range(self.ntorsionmoves):
                # choose bond to rotate around, index is first bead that changes
                index = choose_bond(ca.nrigid - 1, self.P_mid) + 1

                # determine backbone beads
                a1 = np.dot(rotations.aa2mx(ca.rotRigid[index - 1]), np.array([1., 0., 0.]))
                a2 = np.dot(rotations.aa2mx(ca.rotRigid[index]), np.array([1., 0., 0.]))
                x1 = ca.posRigid[index - 1] - 0.4 * a1  # backbone bead
                x2 = ca.posRigid[index] - 0.4 * a2  # backbone bead

                # get bond vector as axis of rotation + random magnitude
                p = x2 - x1
                p /= np.linalg.norm(p)
                p *= 2. * (np.random.random() - 0.5) * self.rotate_backbone
                # convert random rotation to a matrix
                mx = rotations.aa2mx(p)
                # center of rotation is in middle of backbone bond
                center = 0.5 * (x1 + x2)

                # apply rotation to positions and orientations
                for i in range(index, ca.nrigid):
                    a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
                    ca.rotRigid[i] = rotations.rotate_aa(ca.rotRigid[i], p)
                    x = ca.posRigid[i] - 0.4 * a
                    ca.posRigid[i] = np.dot(mx, x - center) + center
                    a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
                    ca.posRigid[i] += 0.4 * a

                # this is necessary for adaptive step taking

    def scale(self, factor):
        self.rotate_backbone *= factor
        self.rotate_base *= factor

    @property
    def stepsize(self):
        return [self.rotate_backbone, self.rotate_base]


# this class should generate a fully random configuration
class OXDNAReseed(takestep.TakestepInterface):
    def __init__(self, radius=3.0):
        self.radius = radius

    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size / 6, coords=coords)

        # random displacement for positions
        ca.posRigid[:] = 2. * self.radius * (np.random.random(ca.posRigid.shape) - 0.5)

        # random rotation for angle-axis vectors
        for rot in ca.rotRigid:
            rot[:] = rotations.random_aa()


# this class should generate a fully random configuration
class OXDNAReseedRandomwalk(takestep.TakestepInterface):
    def __init__(self, radius=3.0):
        self.radius = radius

    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size / 6, coords=coords)

        backbone = np.zeros(3)

        # random rotation for angle-axis vectors
        for pos, rot in zip(ca.posRigid, ca.rotRigid):
            backbone += rotations.vec_random() * 0.7525

            # choose a random rotation
            rot[:] = rotations.random_aa()

            # calcualte center of base from backgone
            a1 = np.dot(rotations.aa2mx(rot), np.array([1., 0., 0.]))
            pos[:] = backbone + 0.4 * a1


def export_xyz(fl, coords):
    ca = CoordsAdapter(nrigid=coords.size / 6, coords=coords)
    fl.write("%d\n\n" % (2 * ca.nrigid))
    for i in range(ca.nrigid):
        a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
        x_back = ca.posRigid[i] - 0.4 * a  # backbone bead
        x_stack = ca.posRigid[i] + 0.4 * a

        fl.write("C %f %f %f\n" % (x_back[0], x_back[1], x_back[2]))
        fl.write("H %f %f %f\n" % (x_stack[0], x_stack[1], x_stack[2]))
