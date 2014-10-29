#include "pele/rotations.h"
#include "pele/aatopology.h"
#include "pele/vecn.h"

namespace pele{

using pele::norm;
using pele::dot;

pele::Array<double>
pele::RigidFragment::to_atomistic(pele::Array<double> const com,
        pele::VecN<3> const & p)
{
    assert(com.size() == _ndim);
    assert(p.size() == 3);
    auto rmat = pele::aa_to_rot_mat(p);
    Array<double> pos(_atom_positions.size());
    HackyMatrix<double> mpos(pos, _ndim);

    // in python this is:
    //      return com + np.dot(R, np.transpose(self.atom_positions)).transpose()
    for (size_t atom = 0; atom<_natoms; ++atom) {
        for (size_t j = 0; j<_ndim; ++j) {
            double val = com[j];
            for (size_t k = 0; k<_ndim; ++k) {
                val += rmat(j,k) * _atom_positions_matrix(atom,k);
            }
            mpos(atom, j) = val;
        }
    }
    return pos;
}

void
pele::RigidFragment::transform_grad(
        pele::VecN<3> const & p,
        pele::Array<double> const g,
        pele::VecN<3> & g_com,
        pele::VecN<3> & g_rot
        )
{
    assert(g.size() == natoms() * 3);
    // view the array as a matrix
    HackyMatrix<double> gmat(g, 3);

    // compute the rotation matrix and derivatives
    pele::MatrixNM<3,3> rmat;
    pele::MatrixNM<3,3> drm1;
    pele::MatrixNM<3,3> drm2;
    pele::MatrixNM<3,3> drm3;
    rot_mat_derivatives(p, rmat, drm1, drm2, drm3);

    // do the center of mass coordinates
    for (size_t k=0; k<3; ++k) {
        double val = 0;
        for (size_t atom=0; atom < _natoms; ++atom) {
            val += gmat(atom,k);
        }
        g_com[k] = val;
    }

    // now do the rotations
    g_rot.assign(0);
    for (size_t atom=0; atom < _natoms; ++atom) {
        double val1 = 0;
        double val2 = 0;
        double val3 = 0;
        for (size_t i=0; i<3; ++i) {
            for (size_t j=0; j<3; ++j) {
                val1 += gmat(atom,i) * drm1(i,j) * _atom_positions_matrix(atom,j);
                val2 += gmat(atom,i) * drm2(i,j) * _atom_positions_matrix(atom,j);
                val3 += gmat(atom,i) * drm3(i,j) * _atom_positions_matrix(atom,j);
            }
        }
        g_rot[0] += val1;
        g_rot[1] += val2;
        g_rot[2] += val3;
    }
}

double
pele::RigidFragment::distance_squared(pele::VecN<3> const & com1, pele::VecN<3> const & p1,
        pele::VecN<3> const & com2, pele::VecN<3> const & p2) const
{
    VecN<3> drij = get_smallest_rij(com1, com2);
    pele::MatrixNM<3,3> R1 = pele::aa_to_rot_mat(p1);
    pele::MatrixNM<3,3> R2 = pele::aa_to_rot_mat(p2);

    MatrixNM<3,3> dR = R2 - R1;  

    double d_M = m_W * dot(drij, drij);
    // we only need the trace, so this can be sped up
    double d_P = dot<3,3,3>(dR, dot<3,3,3>(m_S, transpose<3>(dR))).trace();
    double d_mix = 2. * m_W * dot<3>(drij, dot<3,3>(dR, m_cog));

    double dist2 = d_M + d_P + d_mix;
    return dist2;
}

// sn402 version: overloaded to include a boxlength vector
double
pele::RigidFragment::distance_squared_bulk(pele::VecN<3> const & com1, pele::VecN<3> const & p1,
        pele::VecN<3> const & com2, pele::VecN<3> const & p2, pele::VecN<3> const & boxvec) const
{
    VecN<3> drij = get_smallest_rij_bulk(com1, com2, boxvec);
    pele::MatrixNM<3,3> R1 = pele::aa_to_rot_mat(p1);
    pele::MatrixNM<3,3> R2 = pele::aa_to_rot_mat(p2);

    MatrixNM<3,3> dR = R2 - R1;

    double d_M = m_W * dot(drij, drij);
    // we only need the trace, so this can be sped up
    double d_P = dot<3,3,3>(dR, dot<3,3,3>(m_S, transpose<3>(dR))).trace();
    double d_mix = 2. * m_W * dot<3>(drij, dot<3,3>(dR, m_cog));

    double dist2 = d_M + d_P + d_mix;
    return dist2;
}


void
pele::RigidFragment::distance_squared_grad(pele::VecN<3> const & com1, pele::VecN<3> const & p1,
        pele::VecN<3> const & com2, pele::VecN<3> const & p2,
        VecN<3> & g_M, VecN<3> & g_P
        ) const
{
    VecN<3> drij = get_smallest_rij(com1, com2);
    auto R2 = pele::aa_to_rot_mat(p2);
    MatrixNM<3,3> R1, R11, R12, R13;
    pele::rot_mat_derivatives(p1, R1, R11, R12, R13);

    auto dR = R2 - R1;

    g_M = drij;
    g_M *= -2. * m_W;

    // this linear algebra can be done more efficiently
    auto dRT = pele::transpose(dR);
    g_P[0] = -2. * dot<3,3,3>(R11, dot<3,3,3>(m_S, dRT)).trace();
    g_P[1] = -2. * dot<3,3,3>(R12, dot<3,3,3>(m_S, dRT)).trace();
    g_P[2] = -2. * dot<3,3,3>(R13, dot<3,3,3>(m_S, dRT)).trace();

    // this can also be done more efficiently
    auto temp = dot<3,3>(dR, m_cog);
    temp *= 2. * m_W;
    g_M -= temp;
    g_P[0] -= 2. * m_W * dot<3>(drij, dot<3,3>(R11, m_cog));
    g_P[1] -= 2. * m_W * dot<3>(drij, dot<3,3>(R12, m_cog));
    g_P[2] -= 2. * m_W * dot<3>(drij, dot<3,3>(R13, m_cog));
}

// sn402 version: overloaded to include a boxlength vector
void
pele::RigidFragment::distance_squared_grad_bulk(pele::VecN<3> const & com1, pele::VecN<3> const & p1,
        pele::VecN<3> const & com2, pele::VecN<3> const & p2,
        VecN<3> & g_M, VecN<3> & g_P, pele::VecN<3> const & boxvec
        ) const
{
    VecN<3> drij = get_smallest_rij_bulk(com1, com2, boxvec);
    auto R2 = pele::aa_to_rot_mat(p2);
    MatrixNM<3,3> R1, R11, R12, R13;
    pele::rot_mat_derivatives(p1, R1, R11, R12, R13);

    auto dR = R2 - R1;

    g_M = drij;
    g_M *= -2. * m_W;

    // this linear algebra can be done more efficiently
    auto dRT = pele::transpose(dR);
    g_P[0] = -2. * dot<3,3,3>(R11, dot<3,3,3>(m_S, dRT)).trace();
    g_P[1] = -2. * dot<3,3,3>(R12, dot<3,3,3>(m_S, dRT)).trace();
    g_P[2] = -2. * dot<3,3,3>(R13, dot<3,3,3>(m_S, dRT)).trace();

    // this can also be done more efficiently
    auto temp = dot<3,3>(dR, m_cog);
    temp *= 2. * m_W;
    g_M -= temp;
    g_P[0] -= 2. * m_W * dot<3>(drij, dot<3,3>(R11, m_cog));
    g_P[1] -= 2. * m_W * dot<3>(drij, dot<3,3>(R12, m_cog));
    g_P[2] -= 2. * m_W * dot<3>(drij, dot<3,3>(R13, m_cog));
}


void pele::MeasureAngleAxisCluster::
align(pele::Array<double> const x1, pele::Array<double> x2)
{
    auto c1 = m_topology->get_coords_adaptor(x1);
    auto c2 = m_topology->get_coords_adaptor(x2);

    // now account for the symmetries
    for (size_t isite = 0; isite < m_topology->nrigid(); ++isite) {
        auto const & rotations = m_topology->get_sites()[isite].get_symmetry_rotations();
        auto p1 = c1.get_rb_rotation(isite);
        auto p2 = c2.get_rb_rotation(isite);
        auto mx2 = pele::aa_to_rot_mat(p2);
        auto mx1 = pele::transpose(pele::aa_to_rot_mat(p1));
        auto mx =  pele::dot(mx1, mx2);
        double theta_min = 10.;
        MatrixNM<3,3> rot_best = pele::identity<3>();
        for (auto const & rot : rotations){
            auto mx_diff = dot(mx, rot);
            double theta = norm<3>(rot_mat_to_aa(mx_diff));
            theta -= int(theta / (2. * M_PI)) * 2. * M_PI;
            if (theta < theta_min) {
                theta_min = theta;
                rot_best = rot;
            }
        }
        auto newp2 = rotate_aa(rot_mat_to_aa(rot_best), p2);
        std::copy(newp2.begin(), newp2.end(), p2.begin());
    }
}


Array<double>
pele::RBTopology::to_atomistic(Array<double> rbcoords)
{
    if ( rbcoords.size() != nrigid() * 6 ) {
        throw std::invalid_argument("rbcoords has the wrong size");
    }

    size_t const nrigid = _sites.size();
    auto ca = get_coords_adaptor(rbcoords);
    Array<double> atomistic(3 * natoms_total());
    // view the atomistic coords as a matrix
    size_t istart = 0;
    for (size_t isite=0; isite<nrigid; ++isite) {
        VecN<3> psite = ca.get_rb_rotation(isite);
        auto site_atom_positions = _sites[isite].to_atomistic(
                ca.get_rb_position(isite),
                psite);
        Array<double> atomistic_view(atomistic.view(istart, istart + site_atom_positions.size()));
        atomistic_view.assign(site_atom_positions);

        istart += site_atom_positions.size();
    }
    assert(istart == natoms_total() * 3);
    return atomistic;
}

void
pele::RBTopology::transform_gradient(pele::Array<double> rbcoords,
        pele::Array<double> grad, pele::Array<double> rbgrad)
{
    if ( rbcoords.size() != nrigid() * 6 ) {
        throw std::invalid_argument("rbcoords has the wrong size");
    }
    if (grad.size() != natoms_total() * 3) {
        throw std::invalid_argument("grad has the wrong size");
    }
    if (rbgrad.size() != rbcoords.size()) {
        throw std::invalid_argument("rbgrad has the wrong size");
    }

    CoordsAdaptor ca(nrigid(), 0, rbcoords);
    pele::Array<double> coords_rot(ca.get_rb_rotations());
//        pele::Array<double> rbgrad(rbcoords.size());
    CoordsAdaptor rbgrad_ca(nrigid(), 0, rbgrad);

    size_t istart = 0;
    for (size_t isite=0; isite<nrigid(); ++isite) {
        size_t const site_ndof = _sites[isite].natoms() * 3;
//            std::cout << grad.size() << " " << istart << " " << site_ndof << " " << istart + site_ndof << "\n";
        Array<double> g_site     = grad.view(istart, istart + site_ndof);
        Array<double> p          = ca.get_rb_rotation(isite);
        Array<double> g_com_site = rbgrad_ca.get_rb_position(isite);
        Array<double> g_rot_site = rbgrad_ca.get_rb_rotation(isite);
        _sites[isite].transform_grad(
                ca.get_rb_rotation(isite),
                grad.view(istart, istart + site_ndof),
                rbgrad_ca.get_rb_position(isite),
                rbgrad_ca.get_rb_rotation(isite));
//                p, g_site, g_com_site, g_rot_site);
        istart += site_ndof;
    }
}

pele::VecN<3>
pele::RBTopology::align_angle_axis_vectors(pele::VecN<3> const & p1,
        pele::VecN<3> const & p2in)
{
    pele::VecN<3> p2 = p2in;
    pele::VecN<3> n2, p2n;
    if (norm<3>(p2) < 1e-6) {
        if (norm<3>(p1) < 1e-6) {
            return p2;
        }
        n2 = p1;
        n2 *= 2. * M_PI / norm<3>(p1);
    } else {
        n2 = p2;
        n2 *= 2. * M_PI / norm<3>(p2);
    }

    while (true) {
        p2n = p2;
        p2n += n2;
        if (norm<3>(p2n - p1) > norm<3>(p2 - p1)) {
            break;
        }
        p2 = p2n;
    }

    while (true) {
        p2n = p2;
        p2n -= n2;
        if (norm<3>(p2n - p1) > norm<3>(p2 - p1)) {
            break;
        }
        p2 = p2n;
    }
    return p2;
}

void
pele::RBTopology::align_all_angle_axis_vectors(pele::Array<double> x1,
        pele::Array<double> x2)
{
    auto c1 = get_coords_adaptor(x1);
    auto c2 = get_coords_adaptor(x2);
    for (size_t isite = 0; isite < nrigid(); ++isite) {
        VecN<3> p1 = c1.get_rb_rotation(isite);
        pele::Array<double> p2 = c2.get_rb_rotation(isite);
        auto p2new = align_angle_axis_vectors(p1, p2);
        std::copy(p2new.begin(), p2new.end(), p2.begin());
    }
}

void
pele::RBTopology::align_path(std::list<pele::Array<double> > path)
{
    auto iter1 = path.begin();
    auto iter2 = path.begin();
    iter2++;
    while (iter2 != path.end()) {
        align_all_angle_axis_vectors(*iter1, *iter2);
        ++iter1;
        ++iter2;
    }
}

void
pele::TransformAACluster::rotate(pele::Array<double> x,
        pele::MatrixNM<3,3> const & mx)
{
    auto ca = m_topology->get_coords_adaptor(x);
    if(m_topology->nrigid() > 0) {
        // rotate the center of mass positions by mx
        pele::HackyMatrix<double> rb_pos(ca.get_rb_positions(), 3);
        // make a HackyMatrix view of the transposed rotation matrix
        auto mxT = pele::transpose(mx);
        pele::HackyMatrix<double> mxT_view(mxT.data(), 3, 3);
        assert(mxT_view(0,1) == mx(1,0));
        // do the multiplication
        auto result = hacky_mat_mul(rb_pos, mxT_view);
        // copy the results back into the coordinates array
//        std::cout << "result " << result << std::endl;
        rb_pos.assign(result);

        // rotate each aa rotation by mx
        VecN<3> dp = pele::rot_mat_to_aa(mx);
        for (size_t isite = 0; isite < m_topology->nrigid(); ++isite) {
            pele::Array<double> pview = ca.get_rb_rotation(isite);
            VecN<3> p = pele::rotate_aa(pview, dp);
            // copy the vector back into pview
            std::copy(p.begin(), p.end(), pview.begin());
        }
    }
    if (m_topology->number_of_non_rigid_atoms() > 0) {
        throw std::runtime_error("non-rigid atoms is not yet supported");
//            ca.posAtom[:] = np.dot(mx, ca.posAtom.transpose()).transpose()
    }
}

} // namespace pele
