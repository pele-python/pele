#ifndef _NGT_HPP_
#define _NGT_HPP_
/*
 *
 * This implements the New Graph Transformation method (NGT) described in
 *
 * David Wales, J. Chem. Phys., 2009 http://dx.doi.org/10.1063/1.3133782
 *
 * This procedure computes transition rates and committor probabilities for
 * transition network (kinetic monte carlo).
 */


#include <cstdlib>
#include <iostream>
#include <list>
#include <queue>
#include <assert.h>
#include <stdexcept>
#include <memory>

#include "graph.hpp"

using std::cout;

namespace pele
{

bool compare_degree(node_ptr u, node_ptr v){
    return u->in_out_degree() < v->in_out_degree();
}

class NGT {
public:
    typedef std::map<std::pair<node_id, node_id>, double> rate_map_t;

    std::shared_ptr<Graph> _graph;
    std::set<node_ptr> _A; // the source nodes
    std::set<node_ptr> _B; // the sink nodes
    std::list<node_ptr> intermediates; //this will an up to date list of nodes sorted by the node degree
    bool debug;

    /**
     * the initial waiting time before any graph transformation.  Used to
     * compute steady state rates.
     */
    std::map<node_id, double> initial_tau;
    /**
     * Final values of 1-Pxx for node x after the graph transformation.
     */
    std::map<node_id, double> final_omPxx;
    /**
     * Final values of tau for node x after the graph transformation.
     */
    std::map<node_id, double> final_tau;
    std::map<node_id, double> final_committors;
    std::map<node_id, double> weights; // normally these are equilibrium occupation probabilities


    
    ~NGT()
    {
    }

    /*
     * construct the NGT from an existing graph.
     *
     * The graph will be used directly, without copying.  Any modifications
     * will be reflected in the passed graph
     */
    template<class Acontainer, class Bcontainer>
    NGT(std::shared_ptr<Graph> graph, Acontainer const &A, Bcontainer const &B) :
        _graph(graph),
        debug(false)
    {
        for (auto u : A){
            _A.insert(_graph->get_node(u));
        }
        for (auto u : B){
            _B.insert(_graph->get_node(u));
        }

        // make intermediates
        for (auto const & mapval : _graph->node_map_){
            node_ptr u = mapval.second;
            if (_A.find(u) == _A.end() and _B.find(u) == _B.end()){
                intermediates.push_back(u);
            }
        }

//        std::cout << "number of nodes " << _graph->number_of_nodes() << "\n";
//        std::cout << "A.size() " << _A.size() << "\n";
//        std::cout << "B.size() " << _B.size() << "\n";
//        std::cout << "intermediates.size() " << intermediates.size() << "\n";
        assert(intermediates.size() + _A.size() + _B.size() == _graph->number_of_nodes());

    }

    void set_debug() { debug=true; }
    std::map<node_id, double> const & get_committors() { return final_committors; }

    /*
     * construct the NGT from a map of rate constants.
     */
    template<class Acontainer, class Bcontainer>
    NGT(rate_map_t &rate_constants, Acontainer const &A, Bcontainer const &B) :
        _graph(new Graph()),
        debug(false)
    {
        std::set<node_ptr> nodes;

        // add nodes to the graph and sum the rate constants for all out edges for each node.
        std::map<node_ptr, double> sum_out_rates;
        for (auto const & mapvals : rate_constants){
            node_ptr u = _graph->add_node(mapvals.first.first);
            node_ptr v = _graph->add_node(mapvals.first.second);
            double k = mapvals.second;
            nodes.insert(u);
            nodes.insert(v);

            try {
                sum_out_rates.at(u) += k;
            } catch (std::out_of_range & e) {
                sum_out_rates[u] = k;
            }
        }

        // set tau_x for each node
        // add edge Pxx for each node and initialize P to 0.
        for (auto x : nodes){
            double tau_x = 1. / sum_out_rates[x];
            set_tau(x, tau_x);
            initial_tau[x->id()] = tau_x;
            edge_ptr xx = _graph->_add_edge(x, x);
            set_P(xx, 0.);
        }

        // set Puv for each edge
        for (auto const & mapval : rate_constants){
            node_ptr u = _graph->get_node(mapval.first.first);
            node_ptr v = _graph->get_node(mapval.first.second);
            double k = mapval.second;

            edge_ptr uv = _graph->_add_edge(u, v);
            double tau_u = get_tau(u);
            double Puv = k * tau_u;
            set_P(uv, Puv);

            try {
                sum_out_rates.at(u) += k;
            } catch (std::out_of_range & e) {
                sum_out_rates[u] = k;
            }
        }


        // make the set of A and B
        for (auto a : A){
            _A.insert(_graph->get_node(a));
        }
        for (auto b : B){
            _B.insert(_graph->get_node(b));
        }

        // make a list of intermediates
        for (auto a : _A){
            nodes.erase(a);
        }
        for (auto b : _B){
            nodes.erase(b);
        }
        intermediates.assign(nodes.begin(), nodes.end());


//        std::cout << _graph->number_of_nodes() << "\n";
//        std::cout << _A.size() << "\n";
//        std::cout << _B.size() << "\n";
//        std::cout << intermediates.size() << "\n";
//        std::cout << nodes.size() << "\n";
        assert(intermediates.size() + _A.size() + _B.size() == _graph->number_of_nodes());
    }

    void set_node_occupation_probabilities(std::map<node_id, double> &Peq){
        weights.insert(Peq.begin(), Peq.end());
    }

    /*
     * Sort the list of intermediates.
     *
     * This is done because it is faster to remove nodes with fewer connections first.
     */
    void sort_intermediates(){
        node_ptr x = *intermediates.begin();
        if (debug){
            std::cout << "smallest node degree " << x->in_out_degree() << "\n";
        }
        if (x->in_out_degree() > 4) {
            intermediates.sort(compare_degree);
        }
    }
    
    /*
     * accessors for graph properties P and tau attached to the edges and nodes.
     */
    inline double get_tau(node_ptr u){ return u->tau; }
    inline double get_P(edge_ptr edge){ return edge->P; }
    inline void set_tau(node_ptr u, double tau){ u->tau = tau; }
    inline void set_P(edge_ptr edge, double P){ edge->P = P; }

    /*
     * This returns P for the edge u->u.  This is slow because the edge must first be found.
     */
    double get_node_P(node_ptr u){ return get_P(u->get_successor_edge(u)); }

    /*
     * This returns 1.-P for the edge u->u.
     *
     * If P is close to one compute 1.-P directly by summing P over all the out edges of u.
     * This is extremely important for numerical precision.  It is this ability to deal precisely
     * with both P and 1.-P that makes this method more stable then linear algebra methods.
     */
    double get_node_one_minus_P(node_ptr u){
        edge_ptr uu = u->get_successor_edge(u);
        double Puu = get_P(uu);
        if (Puu < 0.99){
            return 1. - Puu;
        } else {
            // sum the contributions from all other edges
            double omPuu = 0.;
            for (auto eiter = u->out_edge_begin(); eiter != u->out_edge_end(); ++eiter){
                node_ptr v = (*eiter)->head();
                if (v != u){
                    omPuu += (*eiter)->P;
                }
            }
            return omPuu;
        }
    }


    /*
     * node x is being deleted, so update tau for node u
     *
     * tau_u -> tau_u + Pux * tau_x / (1-Pxx)
     */
    void update_node(edge_ptr ux, double omPxx, double tau_x){
        node_ptr u = ux->tail();
        double Pux = get_P(ux);
        double tau_u = get_tau(u);
        double new_tau_u = tau_u + Pux * tau_x / omPxx;
        if (debug){
            std::cout << "updating node " << u->id() << " tau " << tau_u << " -> " << new_tau_u << "\n";
        }
        set_tau(u, new_tau_u);
    }

    /*
     * add an edge to the graph and set P to 0
     */
    edge_ptr add_edge(node_ptr u, node_ptr v){
       edge_ptr edge = _graph->_add_edge(u, v);
       set_P(edge, 0.);
       return edge;
    }

    /*
     * Node x is being deleted, so update P for the edge u -> v
     *
     * Puv -> Puv + Pux * Pxv / (1-Pxx)
     */
    void update_edge(node_ptr u, node_ptr v, edge_ptr ux, edge_ptr xv, double omPxx){
        edge_ptr uv = u->get_successor_edge(v);  // this is slow
        if (uv == NULL){
            uv = add_edge(u, v);
        }

        double Pux = get_P(ux);
        double Pxv = get_P(xv);
        double Puv = get_P(uv);

        double newPuv = Puv + Pux * Pxv / omPxx;
        if (debug) {
            std::cout << "updating edge " << u->id() << " -> " << v->id() << " Puv " << Puv << " -> " << newPuv
                    << " 1-Pxx " << omPxx
                    << " Pux " << Pux
                    << " Pxv " << Pxv
                    << "\n";
        }
        set_P(uv, newPuv);
    }

    /*
     * remove node x from the graph and update its neighbors
     */
    void remove_node(node_ptr x){
        if (debug){
            std::cout << "removing node " << x->id() << "\n";
        }
        double taux = get_tau(x);
//        double Pxx = get_node_P(x);
        double omPxx = get_node_one_minus_P(x);

        // update the node data for all the neighbors
        for (auto eiter = x->in_edge_begin(); eiter != x->in_edge_end(); eiter++){
            edge_ptr edge = *eiter;
            if (edge->tail() != edge->head()){
                update_node(edge, omPxx, taux);
            }
        }

        std::set<node_ptr> neibs = x->in_out_neighbors();
        neibs.erase(x);

        //
        for (auto uxiter = x->in_edge_begin(); uxiter != x->in_edge_end(); ++uxiter){
            edge_ptr ux = *uxiter;
            node_ptr u = ux->tail();
            if (u == x) continue;
            for (auto xviter = x->out_edge_begin(); xviter != x->out_edge_end(); ++xviter){
                edge_ptr xv = *xviter;
                node_ptr v = xv->head();
                if (v == x) continue;
//                if (u == v){
//                    continue;
//                }
                update_edge(u, v, ux, xv, omPxx);
            }
        }

        // remove the node from the graph
        _graph->_remove_node(x);

    }

    /*
     * remove all intermediates from the graph
     */
    void remove_intermediates(){
        while (intermediates.size() > 0){
            sort_intermediates();

            node_ptr x = intermediates.front();
            intermediates.pop_front();

            remove_node(x);
        }
    }

    /*
     * phase one of the rate calculation is to remove all intermediate nodes
     */
    void phase_one(){
        remove_intermediates();
    }

    /*
     * Compute final_tau and final_omPxx for each node in to_remove
     *
     * For each node x in to_remove, this involves removing all other nodes in to_remove, and
     * getting the results from this reduced graph.
     */
    void reduce_all_in_group(std::set<node_ptr> &to_remove, std::set<node_ptr> & to_keep){
        std::list<node_id> Aids, Bids;
        // copy the ids of the nodes in to_remove into Aids
        for (auto u : to_remove){
            Aids.push_back(u->id());
        }
        // copy the ids of the nodes in to_keep into Bids
        for (auto u : to_keep){
            Bids.push_back(u->id());
        }

        // note: should we sort the minima in to_remove?

        if (Aids.size() > 1){
            // make a copy of _graph called working_graph
            auto working_graph = std::make_shared<Graph> (*_graph);
            std::list<node_id> empty_list;
            // make an ngt object for working_graph
            NGT working_ngt(working_graph, std::list<node_id>(), Bids);
            while (Aids.size() > 1){
                /*
                 * Create a new graph and a new NGT object new_ngt.  Pass x as A and Bids as B.  new_ngt will
                 * remove all `intermediates`, i.e. everything in Aids except x.  Then save the final
                 * value of 1-Pxx and tau_x.
                 */
                // choose an element x and remove it from the list
                node_id x = Aids.back();
                Aids.pop_back();
                std::list<node_id> newAids;
                newAids.push_back(x);

                // make a new graph from the old graph
                auto new_graph = std::make_shared<Graph>(*working_graph);

                // remove all nodes from new_graph except x
                NGT new_ngt(new_graph, newAids, Bids);
                new_ngt.remove_intermediates();
                node_ptr xptr = new_graph->get_node(x);
                final_omPxx[x] = new_ngt.get_node_one_minus_P(xptr);
                final_tau[x] = new_ngt.get_tau(xptr);

                // delete node x from the old_graph
                working_ngt.remove_node(working_graph->get_node(x));
            }
            // there is one node left. we can just read off the results
            assert(Aids.size() == 1);
            node_id x = Aids.back();
            Aids.pop_back();
            node_ptr xptr = working_graph->get_node(x);
            final_omPxx[x] = working_ngt.get_node_one_minus_P(xptr);
            final_tau[x] = working_ngt.get_tau(xptr);

        } else if (Aids.size() == 1) {
            // if there is only one node in A then we can just read off the results.
            node_id x = Aids.back();
            Aids.pop_back();
            node_ptr xptr = _graph->get_node(x);
            final_omPxx[x] = get_node_one_minus_P(xptr);
            final_tau[x] = get_tau(xptr);
        }
        assert(Aids.size() == 0);
    }

    /*
     * Phase two, compute final_tau and final_omPxx for each x separately in _A and in _B
     */
    void phase_two(){
        reduce_all_in_group(_A, _B);
        reduce_all_in_group(_B, _A);
    }

    /*
     * do phase one and phase two of the rate calculation
     */
    void compute_rates(){
        phase_one();
        phase_two();
    }

    /*
     * compute the final rate A->B or B->A from final_tau and final_omPxx
     */
    double _get_rate_final(std::set<node_ptr> &A){
        double rate_sum = 0.;
        double norm = 0.;
        for (auto a : A){
            double omPxx = final_omPxx.at(a->id());
            double tau_a = final_tau.at(a->id());
            double weight = 1.;
            if (weights.size() > 0){
                weight = weights.at(a->id());
            }
            rate_sum += weight * omPxx / tau_a;
            norm += weight;
        }
        return rate_sum / norm;
    }

    /*
     * Return the rate A->B
     */
    double get_rate_AB(){
        return _get_rate_final(_A);
    }

    /*
     * Return the rate B->A
     */
    double get_rate_BA(){
        return _get_rate_final(_B);
    }

    double _get_rate_SS(std::set<node_ptr> & A, std::set<node_ptr> & B){
        double kAB = 0.;
        double norm = 0.;
        for (auto a : A){
            // compute PaB the probability that this node goes directly to B
            double PaB = 0.;
            for (auto eiter = a->out_edge_begin(); eiter != a->out_edge_end(); ++eiter){
                edge_ptr ab = *eiter;
                node_ptr b = ab->head();
                if (B.find(b) != B.end()){
                    PaB += get_P(ab);
                }
            }
            double weight = 1.;
            if (weights.size() > 0){
                weight = weights.at(a->id());
            }
            kAB += weight * PaB / initial_tau.at(a->id());
            norm += weight;
        }
        return kAB / norm;
    }

    /*
     * Return the steady state rate A->B
     *
     * this must be called after calling phase_one
     */
    double get_rate_AB_SS(){
        return _get_rate_SS(_A, _B);
    }

    /*
     * Return the steady state rate B->A
     *
     * this must be called after calling phase_one
     */
    double get_rate_BA_SS(){
        return _get_rate_SS(_B, _A);
    }

    /*
     * sum the probabilities of the out edges of x that end in B normalized by 1-Pxx
     */
    double get_PxB(node_ptr x, std::set<node_id> & B){
        double PxB = 0.;
        double Pxx = 0.;
        double omPxx = 0.;
        for (auto eiter = x->out_edge_begin(); eiter != x->out_edge_end(); ++eiter){
            edge_ptr xb = *eiter;
            node_ptr b = xb->head();
            double Pxb = get_P(xb);
            if (b == x){
                Pxx = Pxb;
            } else {
                omPxx += Pxb;
            }
            if (B.find(b->id()) != B.end()){
                PxB += Pxb;
            }
        }
        if (Pxx < 0.9){
            omPxx = 1. - Pxx;
        }
        return PxB / omPxx;
    }

    /*
     * compute the committors for all intermediates
     *
     * \param to_remove a list of nodes that will be removed.  Committor values will
     *     be computed for these nodes
     * \param to_keep a list of nodes that should not be deleted.
     * \param committor_targets a list of nodes that should not be deleted.  These
     *     nodes will be the targets in the committor calucation.
     *
     * All nodes should be in one of the three passed groups of nodes.  Duplicates
     * between to_keep and committor_targets are OK.
     */
    void _remove_nodes_and_compute_committors(std::list<node_ptr> &to_remove,
            std::set<node_ptr> &to_keep, std::set<node_ptr> const &committor_targets)
    {
        // make a copy of to_remove.  Store the id's
        std::list<node_id> to_remove_cp;
        for (auto u : to_remove){
            to_remove_cp.push_back(u->id());
        }

        // copy the nodes from to_keep and committor_target into a new set Bids;
        // make a copy of committor_target
        std::set<node_id> Bids;
        std::set<node_id> targets;
        // create a set of Bids
        for (auto u : to_keep){
            Bids.insert(u->id());
        }
        for (auto u : committor_targets){
            Bids.insert(u->id());
            targets.insert(u->id());
        }

        // ensure there are no unaccounted for nodes
        assert(to_remove_cp.size() + Bids.size() == _graph->number_of_nodes());

        // note: should we sort the nodes in to_remove?

        while (to_remove_cp.size() > 0){
            /*
             * Create a new graph and a new NGT object new_ngt.  Pass x as A and Bids as B.  new_ngt will
             * remove all `intermediates`, i.e. everything in to_remove except x.  Then save the final
             * value of 1-Pxx and tau_x.
             */
            // choose an element x and remove it from the list
            node_id x = to_remove_cp.back();
            to_remove_cp.pop_back();
            std::list<node_id> Aids;
            Aids.push_back(x);

            // make a copy of _graph
            auto new_graph = std::make_shared<Graph>(*_graph);

            // remove all to_remove nodes from new_graph except x
            NGT new_ngt(new_graph, Aids, Bids);
            new_ngt.remove_intermediates();
            node_ptr xptr = new_graph->get_node(x);
            final_omPxx[x] = new_ngt.get_node_one_minus_P(xptr);
            final_tau[x] = new_ngt.get_tau(xptr);
            if (! targets.empty()){
                final_committors[x] = new_ngt.get_PxB(xptr, targets);
            }

            // delete node x from _graph
            this->remove_node(_graph->get_node(x));
        }
    }

    /*
     * Compute the rate from A->B and committor probabilities for all intermediates
     *
     * This is much slower than compute_rates.
     * If you don't want committors use that function instead
     */
    void compute_rates_and_committors(){
        _remove_nodes_and_compute_committors(intermediates, _A, _B);
        intermediates.clear();

        phase_two();

        // set the committor for nodes in A to 0
        for (auto a : _A){
            final_committors[a->id()] = 0.;
        }
        // set the committor for nodes in B to 1
        for (auto b : _B){
            final_committors[b->id()] = 1.;
        }
    }


};

}
#endif
