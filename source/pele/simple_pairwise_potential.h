#ifndef PYGMIN_SIMPLE_PAIRWISE_POTENTIAL_H
#define PYGMIN_SIMPLE_PAIRWISE_POTENTIAL_H

#include "base_potential.h"
#include "array.h"
#include "distance.h"

namespace pele
{
    /**
     * Define a base class for potentials with simple pairwise interactions that
     * depend only on magnitude of the atom separation
     *
     * This class loops though atom pairs, computes the distances and get's the
     * value of the energy and gradient from the class pairwise_interaction.
     * pairwise_interaction is a passed parameter and defines the actual
     * potential function.
     */
    template<typename pairwise_interaction, 
                 typename distance_policy = cartesian_distance >
    class SimplePairwisePotential : public BasePotential
    {
    protected:
        pairwise_interaction *_interaction;
        distance_policy *_dist;

        SimplePairwisePotential(pairwise_interaction *interaction,
                distance_policy *dist=NULL) : 
            _interaction(interaction), _dist(dist) 
        {
            if(_dist == NULL) _dist = new distance_policy;
        }

    public:
        virtual ~SimplePairwisePotential() 
        { 
            if (_interaction != NULL) delete _interaction; 
            if (_dist != NULL) delete _dist; 
        }

        virtual double get_energy(Array<double> x);
        virtual double get_energy_gradient(Array<double> x, Array<double> grad);
        virtual void get_hessian(Array<double> x, Array<double> hess);
    };

    template<typename pairwise_interaction, typename distance_policy>
    inline double SimplePairwisePotential<pairwise_interaction,distance_policy>::get_energy_gradient(Array<double> x, Array<double> grad)
    {
        double e=0.;
        double gij, dr[3];
        const size_t n = x.size();
        const size_t natoms = x.size()/3;

        for(size_t i=0; i<n; ++i)
            grad[i] = 0.;

        for(size_t atomi=0; atomi<natoms; ++atomi) {
            int i1 = 3*atomi;
            for(size_t atomj=atomi+1; atomj<natoms; ++atomj) {
                int j1 = 3*atomj;

                _dist->get_rij(dr, &x[i1], &x[j1]);

                double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                e += _interaction->energy_gradient(r2, &gij, atomi, atomj);
                for(size_t k=0; k<3; ++k)
                    grad[i1+k] -= gij * dr[k];
                for(size_t k=0; k<3; ++k)
                    grad[j1+k] += gij * dr[k];
            }
        }
        return e;
    }

    template<typename pairwise_interaction, typename distance_policy>
        inline void SimplePairwisePotential<pairwise_interaction,distance_policy>::get_hessian(Array<double> x, Array<double> hess)
        {
    		double hij, gij, dr[3], r2, e;
            const size_t N = x.size();
            const size_t N2 = N*N;
            const size_t natoms = x.size()/3;
            assert(hess.size() == N2);
            hess.assign(0.);


            //LOOP ONLY OVER THE LOWER TRIANGULAR PART OF THE MATRIX
            //build array of distances and interaction
            Array<double*> arraydr(natoms*natoms);
            for(size_t atomi=0; atomi<natoms; ++atomi) {
				int i1 = 3*atomi;
				for(size_t atomj=0; atomj<=atomi; ++atomj) {
					int j1 = 3*atomj;
					_dist->get_rij(dr, &x[i1], &x[j1]);
					double * new_dr = new double[3];
					std::copy(dr,dr+3,new_dr);
					arraydr[natoms*atomi+atomj] = new_dr;
					arraydr[natoms*atomj+atomi] = new_dr;
				}
            }

            for(size_t atomi=0; atomi<natoms; ++atomi) {
				int i1 = 3*atomi;
				//diagonal blocks
				for (size_t atomj=0;atomj<natoms;++atomj){
					if (atomj != atomi){
						int j1 = 3*atomj;
						double* dr = arraydr[natoms*atomi+atomj];
						r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
						_interaction->hessian(r2, &hij, atomi, atomj);
						e += 0.5 * _interaction->energy_gradient(r2, &gij, atomi, atomj);
						for (size_t k=0;k<3;++k){
							//diagonal block - diagonal terms
							double Hii_diag = (hij+gij)*dr[k]*dr[k]/r2 - gij;
							hess[N*(i1+k)+i1+k] += Hii_diag;
							//off diagonal block - diagonal terms
							double Hij_diag = -Hii_diag;
							hess[N*(i1+k)+j1+k] = Hij_diag;
							for (size_t l = k+1;l<3;++l){
								//diagonal block - off diagonal terms
								double Hii_off = (hij+gij)*dr[k]*dr[l]/r2;
								hess[N*(i1+k)+i1+l] += Hii_off;
								hess[N*(i1+l)+i1+k] += Hii_off;
								//off diagonal block - off diagonal terms
								double Hij_off = -Hii_off;
								hess[N*(i1+k)+j1+l] = Hij_off;
								hess[N*(i1+l)+j1+k] = Hij_off;
							}
						}
					}
				}
            }
            /*//LOOP ONLY OVER THE LOWER TRIANGULAR PART OF THE MATRIX
            //we use the fact the gij return the the -derivative/r
            for(size_t atomi=0; atomi<natoms; ++atomi) {
                int i1 = 3*atomi;
                for(size_t atomj=0; atomj<=atomi; ++atomj) {
                    int j1 = 3*atomj;
                    //diagonal blocks
                    if (atomi == atomj)
                    {
                    	for (size_t a=0;a<natoms;++a){
							if (a != atomi){
								double* dr = arraydr[natoms*atomi+a];
								r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
								_interaction->hessian(r2, &hij, atomi, a);
								e = _interaction->energy_gradient(r2, &gij, atomi, a);
								for (size_t k=0;k<3;++k){
									hess[N*(i1+k)+j1+k] += (hij+gij)*dr[k]*dr[k]/r2 - gij;
									for (size_t l = k+1;l<3;++l){
										double Hii_off = (hij+gij)*dr[k]*dr[l]/r2;
										hess[N*(i1+k)+j1+l] += Hii_off;
										hess[N*(i1+l)+j1+k] += Hii_off;
									}
								}
							}
                    	}
                    }
                    //non diagonal blocks
                    else
                    {
						//diagonal terms
						for (size_t k=0;k<3;++k)
						{
							double* dr = arraydr[natoms*atomi+atomj];
							r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
							_interaction->hessian(r2, &hij, atomi, atomj);
							e = _interaction->energy_gradient(r2, &gij, atomi, atomj);
							double Hij_diag = -(gij+hij)*dr[k]*dr[k]/r2 + gij;
							hess[N*(i1+k)+j1+k] = Hij_diag;
							hess[N*(j1+k)+i1+k] = Hij_diag;
							//non diagonal terms
							for (size_t l = k+1;l<3;++l)
							{
								double Hij_off = -(gij+hij)*dr[k]*dr[l]/r2;
								hess[N*(i1+k)+j1+l] = Hij_off;
								hess[N*(i1+l)+j1+k] = Hij_off;
								hess[N*(j1+k)+i1+l] = Hij_off;
								hess[N*(j1+l)+i1+k] = Hij_off;
							}
						}
					}
                }
            }*/

            /*//print hessian
            for (size_t i=0; i<hess.size(); ++i){
            	std::cout<<hess[i]<<" ";
            	int j = i+1;
            	if (j % (int) (3*sqrt(hess.size())) == 0)
					std::cout<<"\n"<<std::endl;
				else if (j % (int) sqrt(hess.size()) == 0)
					std::cout<<""<<std::endl;
				else if (j % 3 == 0){
					std::cout<<"  ";
            	}
            }*/
        }

    template<typename pairwise_interaction, typename distance_policy>
    inline double SimplePairwisePotential<pairwise_interaction, distance_policy>::get_energy(Array<double> x)
    {
        double e=0.;
        size_t const natoms = x.size()/3;

        for(size_t atomi=0; atomi<natoms; ++atomi) {
            size_t i1 = 3*atomi;
            for(size_t atomj=atomi+1; atomj<natoms; ++atomj) {
                size_t j1 = 3*atomj;
                double dr[3];
                _dist->get_rij(dr, &x[i1], &x[j1]);
                double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                e += _interaction->energy(r2, atomi, atomj);
            }
        }
        return e;
    }
}

#endif
