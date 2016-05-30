#ifndef _PELE_XY_MODEL_ONLINE_H
#define _PELE_XY_MODEL_ONLINE_H

#include <fstream>
#include <iostream>
#include <vector>

#include "pele/base_potential_online.h"

namespace pele {

/**
 * Adjacency list graph implementation according to S. S. Skiena,
 * Algorithm Design Manual 
 */

struct edgenode {
    const unsigned int y;
    edgenode* next;
    edgenode(const unsigned int y_, edgenode*const& next_)
        : y(y_),
          next(next_)
    {}
};

class graph {
    std::vector<edgenode*> edges;
    std::vector<unsigned int> degree;
    unsigned int nvertices_;
    unsigned int nedges_;
    const bool directed_;
public:
    graph(const unsigned int maxv_=1000, const bool directed__=false)
        : directed_(directed__)
    {
        resize(maxv_);
    }
    graph(const unsigned int nvertices_input, pele::Array<size_t> x,
        pele::Array<size_t> y)
        : directed_(false)
    {
        resize(nvertices_input);
        loadvec(nvertices_input, x, y);
    }
    ~graph()
    {
        for (unsigned int i = 0; i < nvertices_; ++i) {
            edgenode* p = edges.at(i);
            while (p != NULL) {
                edgenode* tmp = p;
                p = p->next;
                delete tmp;
            }
        }
    }
    void loadtxt(const char* name)
    {
        /**
         * File format (number of lines == nedges + 1): 
         * First Line: nvertices and nedges
         * Following nedges lines: two numbers:
         * 1. vertex A, 2. vertex B
         */
        std::ifstream inp(name);
        unsigned int nvertices_file;
        unsigned int nedges_file;
        inp >> nvertices_file;
        inp >> nedges_file;
        std::cout << "number of vertices in file: " << nvertices_file << "\n";
        std::cout << "number of edges in file: " << nedges_file << "\n";
        for (unsigned int i = 0; i < nedges_file; ++i) {
            unsigned int x;
            unsigned int y;
            inp >> x;
            inp >> y;
            std::cout << "insert edge (" << x << ", " << y << ")\n";
            insert_edge(x, y, directed());
        }
        nvertices_ = nvertices_file;
    }
    void loadvec(const size_t nvertices_input, pele::Array<size_t> x, pele::Array<size_t> y)
    {
        if (x.size() != y.size()) {
            throw std::runtime_error("graph: illegal input: x.size() != y.size()");
        }
        for (size_t i = 0; i < x.size(); ++i) {
            insert_edge(x[i], y[i], directed());
        }
        if (nvertices_input > maxv()) {
            throw std::runtime_error("graph: illegal input: too many vertices");
        }
        nvertices_ = nvertices_input;
    }
    void insert_edge(const unsigned int x, const unsigned int y, bool no_mirror)
    {
        edgenode* p = new edgenode(y, edges.at(x));
        edges.at(x) = p;
        ++degree.at(x);
        if (!no_mirror) {
            insert_edge(y, x, true);
        }
        else {
            ++nedges_;
        }
    }
    void resize(const unsigned int maxv_)
    {
        edges.assign(maxv_, NULL);
        degree.assign(maxv_, 0);
        nvertices_ = 0;
        nedges_ = 0;
    }
    void print(std::ostream& stm)
    {
        stm << "number of vertices in graph: " << nvertices() << "\n";
        stm << "number of edges in graph: " << nedges() << "\n";
        stm << "adjacency list:\n";
        for (unsigned int i = 0; i < nvertices_; ++i) {
            stm << i << " ==>> ";
            edgenode* p = edges.at(i);
            unsigned int tmp = 0;
            while (p != NULL) {
                if (tmp++) {
                    stm << " --> ";
                }
                stm << p->y;
                p = p->next;
            }
            stm << "\n";
        }
    }
    unsigned int nvertices() const
    {
        return nvertices_;
    }
    unsigned int nedges() const
    {
        return nedges_;
    }
    bool directed() const
    {
        return directed_;
    }
    unsigned int maxv() const
    {
        return edges.size();
    }
    edgenode* get_edges(const unsigned int v) const
    {
        return edges.at(v);
    }
};

/**
 * XY model implementation where the energy is computed by summing each
 * linked list in the adjacency list vector in tun. 
 */
class XYModelOnline : public BasePotentialOnline {
    graph m_topology;
public:
    XYModelOnline(const size_t nr_spins, pele::Array<size_t> head, pele::Array<size_t> tail)
        : BasePotentialOnline(nr_spins),
          m_topology(nr_spins, head, tail)
    {}
    double get_energy(Array<double> x, const size_t batch_number)
    {
        /**
         * Return ith term of the potential energy, i.e., the sum of 
         * -cos(x[i] - x[j]) for all j which are linked to i
         */
        if (x.size() != m_topology.nvertices()) {
            throw std::runtime_error("XYModelOnline: x.size() != nr vertices in graph");
        }
        double energy = 0;
        edgenode* en = m_topology.get_edges(batch_number);
        while (en) {
            energy += -std::cos(x[batch_number] - x[en->y]);
            en = en->next;
        }
        return energy;
    }
    double get_energy_gradient_batch(Array<double> x, const size_t batch_number,
        Array<double> ograd)
    {
        /**
         * Return full potential energy. Compute
         * gradient of ith term.
         */
        if (x.size() != m_topology.nvertices() || ograd.size() != x.size()) {
            throw std::runtime_error("XYModelOnline: illegal input");
        }
        ograd.assign(0);
        edgenode* en = m_topology.get_edges(batch_number);
        while (en) {
            const double tmp = std::sin(x[batch_number] - x[en->y]);
            for (size_t k = 0; k < x.size(); ++k) {
                ograd[k] += ((k == batch_number) - (k == en->y)) * tmp;
            }
            en = en->next;
        }
        return BasePotentialOnline::get_energy(x);
    }
    double get_energy_gradient_gradient2_batch(Array<double>x,
        const size_t batch_number, Array<double> ograd, Array<double> ograd2)
    {
        /**
         * Return full potential energy. Compute
         * 2nd gradient of ith term.
         */
        if (x.size() != m_topology.nvertices() || ograd.size() != x.size() || ograd2.size() != x.size()) {
            throw std::runtime_error("XYModelOnline: illegal input");
        }
        ograd.assign(0);
        ograd2.assign(0);
        edgenode* en = m_topology.get_edges(batch_number);
        while (en) {
            const double tmp = std::sin(x[batch_number] - x[en->y]);
            const double tmp2 = std::cos(x[batch_number] - x[en->y]);
            for (size_t k = 0; k < x.size(); ++k) {
                ograd[k] += ((k == batch_number) - (k == en->y)) * tmp;
                ograd2[k] += ((k == batch_number) + (k == en->y)) * tmp2;
            }
            en = en->next;
        }
        return BasePotentialOnline::get_energy(x);
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_XY_MODEL_ONLINE_H

