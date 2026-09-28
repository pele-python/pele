#include <iostream>
#include "pele/lammps_pele.h"
#include "library.h"
#if defined(LAMMPS_LIB_MPI)
#include <mpi.h>
#endif

#include <stdio.h>
#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "info.h"
#include "input.h"
#include "improper.h"
#include "kspace.h"
#include "lmppython.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#if defined(LMP_PLUGIN)
#include "plugin.h"
#endif
#include "region.h"
#include "timer.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include "pele/base_potential.h"
#include <cstdint>

using namespace std;
using namespace LAMMPS_NS;

namespace pele {

LAMMPSPotential::LAMMPSPotential(uintptr_t handle) : handle(handle) {};

double LAMMPSPotential::get_energy(Array<double> xs)
{
    auto lmp = (LAMMPS *) this->handle;
    auto comm = (Comm *) lmp->comm;
    auto timer = (Timer *) lmp->timer;
    auto modify = (Modify *) lmp->modify;
    auto neighbor = (Neighbor *) lmp->neighbor;
    auto domain = (Domain *) lmp->domain;
    auto atom = (Atom *) lmp->atom;
    auto force = (Force *) lmp->force;
    auto update = (Update *) lmp->update;
    auto output = (Output *) lmp->output;
    class Compute *pe_compute;
    double **x = lmp->atom->x;
    int ii;
    int eflag=1;
    int vflag=1;
    int triclinic = domain->triclinic;
    int neigh_every, neigh_delay, neigh_dist_check;

    if (update->first_update == 0) {
        lammps_command(lmp, "run 0 pre no post no");
    }

    neigh_every = neighbor->every;
    neigh_delay = neighbor->delay;
    neigh_dist_check = neighbor->dist_check;

    if ((neigh_every != 1) || (neigh_delay != 0)) {
        neighbor->every = 1;
        neighbor->delay = 0;
        neighbor->dist_check = 1;
    }

    update->whichflag = 1;

    pe_compute = modify->get_compute_by_id("thermo_pe");
    if (!pe_compute) lmp->error->all(FLERR,"Minimization could not find thermo_pe compute");

    // read coordinates
    for (int i = 0; i < lmp->atom->nlocal; i++) {
        ii = 3*i;
        x[i][0] = xs[ii+0];
        x[i][1] = xs[ii+1];
        x[i][2] = xs[ii+2];
    }

    // check for reneighboring
    // always communicate since minimizer moved atoms

    int nflag = neighbor->decide();

    if (nflag == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(Timer::COMM);
    } else {
        if (modify->n_min_pre_exchange) {
            timer->stamp();
            modify->min_pre_exchange();
            timer->stamp(Timer::MODIFY);
        }
        if (triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        if (domain->box_change) {
            domain->reset_box();
            comm->setup();
            if (neighbor->style) neighbor->setup_bins();
        }
        timer->stamp();
        comm->exchange();
        if (atom->sortfreq > 0 &&
                update->ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(Timer::COMM);
        if (modify->n_min_pre_neighbor) {
            modify->min_pre_neighbor();
            timer->stamp(Timer::MODIFY);
        }
        neighbor->build(1);
        timer->stamp(Timer::NEIGH);
        if (modify->n_min_post_neighbor) {
            modify->min_post_neighbor();
            timer->stamp(Timer::MODIFY);
        }
    }

    memset(&atom->f[0][0],0,3*sizeof(double) * atom->nlocal);

    timer->stamp();

    if (modify->n_min_pre_force) {
        modify->min_pre_force(vflag);
        timer->stamp(Timer::MODIFY);
    }

    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);

    if (atom->molecular != Atom::ATOMIC) {
        if (force->bond) force->bond->compute(eflag,vflag);
        if (force->angle) force->angle->compute(eflag,vflag);
        if (force->dihedral) force->dihedral->compute(eflag,vflag);
        if (force->improper) force->improper->compute(eflag,vflag);
        timer->stamp(Timer::BOND);
    }

    if (force->kspace) {
        force->kspace->setup();
        force->kspace->compute(eflag,vflag);
        timer->stamp(Timer::KSPACE);
    }

    if (modify->n_min_pre_reverse) {
        modify->min_pre_reverse(eflag,vflag);
        timer->stamp(Timer::MODIFY);
    }

    if (force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
    }

    if (modify->n_min_post_force) {
         timer->stamp();
         modify->min_post_force(vflag);
         timer->stamp(Timer::MODIFY);
    }

    neighbor->every = neigh_every;
    neighbor->delay = neigh_delay;
    neighbor->dist_check = neigh_dist_check;

    update->whichflag = 0;

    pe_compute->compute_scalar();
    return pe_compute->scalar;
}

double LAMMPSPotential::get_energy_gradient(Array<double> xs, Array<double> gs)
{
    auto lmp = (LAMMPS *) this->handle;
    auto comm = (Comm *) lmp->comm;
    auto timer = (Timer *) lmp->timer;
    auto modify = (Modify *) lmp->modify;
    auto neighbor = (Neighbor *) lmp->neighbor;
    auto domain = (Domain *) lmp->domain;
    auto atom = (Atom *) lmp->atom;
    auto force = (Force *) lmp->force;
    auto update = (Update *) lmp->update;
    auto output = (Output *) lmp->output;
    class Compute *pe_compute;
    double **x = lmp->atom->x;
    double **f = lmp->atom->f;
    int ii;
    int eflag=1;
    int vflag=1;
    int triclinic = domain->triclinic;
    int neigh_every, neigh_delay, neigh_dist_check;

    if (update->first_update == 0) {
        lammps_command(lmp, "run 0 pre no post no");
    }

    neigh_every = neighbor->every;
    neigh_delay = neighbor->delay;
    neigh_dist_check = neighbor->dist_check;

    if ((neigh_every != 1) || (neigh_delay != 0)) {
        neighbor->every = 1;
        neighbor->delay = 0;
        neighbor->dist_check = 1;
    }

    update->whichflag = 1;

    pe_compute = modify->get_compute_by_id("thermo_pe");
    if (!pe_compute) lmp->error->all(FLERR,"Minimization could not find thermo_pe compute");

    // read coordinates
    for (int i = 0; i < lmp->atom->nlocal; i++) {
        ii = 3*i;
        x[i][0] = xs[ii+0];
        x[i][1] = xs[ii+1];
        x[i][2] = xs[ii+2];
    }

    // check for reneighboring
    // always communicate since minimizer moved atoms

    int nflag = neighbor->decide();

    if (nflag == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(Timer::COMM);
    } else {
        if (modify->n_min_pre_exchange) {
            timer->stamp();
            modify->min_pre_exchange();
            timer->stamp(Timer::MODIFY);
        }
        if (triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        if (domain->box_change) {
            domain->reset_box();
            comm->setup();
            if (neighbor->style) neighbor->setup_bins();
        }
        timer->stamp();
        comm->exchange();
        if (atom->sortfreq > 0 &&
                update->ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(Timer::COMM);
        if (modify->n_min_pre_neighbor) {
            modify->min_pre_neighbor();
            timer->stamp(Timer::MODIFY);
        }
        neighbor->build(1);
        timer->stamp(Timer::NEIGH);
        if (modify->n_min_post_neighbor) {
            modify->min_post_neighbor();
            timer->stamp(Timer::MODIFY);
        }
    }

    memset(&atom->f[0][0],0,3*sizeof(double) * atom->nlocal);

    timer->stamp();

    if (modify->n_min_pre_force) {
        modify->min_pre_force(vflag);
        timer->stamp(Timer::MODIFY);
    }

    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);

    if (atom->molecular != Atom::ATOMIC) {
        if (force->bond) force->bond->compute(eflag,vflag);
        if (force->angle) force->angle->compute(eflag,vflag);
        if (force->dihedral) force->dihedral->compute(eflag,vflag);
        if (force->improper) force->improper->compute(eflag,vflag);
        timer->stamp(Timer::BOND);
    }

    if (force->kspace) {
        force->kspace->setup();
        force->kspace->compute(eflag,vflag);
        timer->stamp(Timer::KSPACE);
    }

    if (modify->n_min_pre_reverse) {
        modify->min_pre_reverse(eflag,vflag);
        timer->stamp(Timer::MODIFY);
    }

    if (force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
    }

    if (modify->n_min_post_force) {
         timer->stamp();
         modify->min_post_force(vflag);
         timer->stamp(Timer::MODIFY);
    }

    // write force
    for (int i = 0; i < lmp->atom->nlocal; i++) {
        ii = 3*i;
        gs[ii+0] = -f[i][0];
        gs[ii+1] = -f[i][1];
        gs[ii+2] = -f[i][2];
    }

    neighbor->every = neigh_every;
    neighbor->delay = neigh_delay;
    neighbor->dist_check = neigh_dist_check;

    update->whichflag = 0;

    pe_compute->compute_scalar();
    return pe_compute->scalar;
}

}