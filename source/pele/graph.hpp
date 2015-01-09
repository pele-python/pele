#ifndef _GRAPH_CPP_H_
#define _GRAPH_CPP_H_
/*
 * This implements a lightweight directed graph.  This was primarily designed for
 * use in the New Graph Transformation (NGT) method.  This graph was optimized
 * for the following actions to be as fast as possible
 *
 *     * removing nodes
 *     * iteration over out edges
 *     * iteration over in edges
 *     * return the edge u->v
 *     * access node property `double `
 *     * access edge property `double P`
 *     * allow for loop edges u->u
 *     * copy graph
 *
 * The requirement for fast removing of nodes means I can't assign each node
 * an index and use std::vector's.  This is at odds with fast access to the edge u->v.
 * The solution here is to use a std::map Node.successor_edge_
 */

#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <assert.h>
#include <stdexcept>
#include <memory>


namespace pele
{
class Edge;
class Node;
typedef size_t node_id;
typedef Node * node_ptr;
typedef Edge * edge_ptr;
typedef int color_type;
color_type color_white = 0;
color_type color_grey = 1;
color_type color_black = 4;

/**
 * basic class for an edge (arc) in the graph
 */
class Edge{
public:
    node_ptr head_; // node the edge points to
    node_ptr tail_; // node the edge comes from
    double P;

    Edge(node_ptr tail, node_ptr head):
        head_(head),
        tail_(tail)
    {}

    node_ptr head(){ return head_; }
    node_ptr tail(){ return tail_; }
};

/**
 * basic class for a node in the graph
 */
class Node{
public:
    typedef std::set<edge_ptr> edge_list;
    typedef typename edge_list::iterator edge_iterator;
    typedef std::map<node_ptr, edge_ptr> successor_map_t;
private:
    node_id id_;
    edge_list out_edge_list_; // list of outgoing edges
    edge_list in_edge_list_; // list of outgoing edges


public:
    successor_map_t successor_map_;
    double tau;

    Node(node_id id):
        id_(id)
    {}

    void add_out_edge(edge_ptr edge){
        out_edge_list_.insert(edge);
        successor_map_[edge->head()] = edge;
    }
    void remove_out_edge(edge_ptr edge){
        out_edge_list_.erase(edge);
        successor_map_.erase(edge->head());
    }
    edge_list & get_out_edges(){ return out_edge_list_; }
    edge_iterator out_edge_begin(){ return out_edge_list_.begin(); }
    edge_iterator out_edge_end(){ return out_edge_list_.end(); }

    void add_in_edge(edge_ptr edge){ in_edge_list_.insert(edge); }
    void remove_in_edge(edge_ptr edge){ in_edge_list_.erase(edge); }
    edge_list & get_in_edges(){ return in_edge_list_; }
    edge_iterator in_edge_begin(){ return in_edge_list_.begin(); }
    edge_iterator in_edge_end(){ return in_edge_list_.end(); }

    node_id id() const { return id_; }
    size_t out_degree() const { return out_edge_list_.size(); }
    size_t in_degree() const { return in_edge_list_.size(); }
    size_t in_out_degree() const { return out_degree() + in_degree(); }
    std::set<node_ptr> in_out_neighbors();

    /*
     * return the edge u->v
     */
    edge_ptr get_successor_edge(node_ptr v){
        successor_map_t::iterator miter = successor_map_.find(v);
        if (miter == successor_map_.end()){
            return NULL;
        } else{
            return miter->second;
        }
    }
};

std::set<node_ptr> Node::in_out_neighbors() {
    std::set<node_ptr> neibs;
    Node::edge_iterator eiter;
    for (eiter = out_edge_begin(); eiter != out_edge_end(); eiter++){
        neibs.insert((*eiter)->head());
    }
    for (eiter = in_edge_begin(); eiter != in_edge_end(); eiter++){
        neibs.insert((*eiter)->tail());
    }
    return neibs;
}


class Graph
{
public:
    typedef std::map<node_id, node_ptr> node_map_t;
    node_map_t node_map_;
    typedef std::set<edge_ptr> edge_list_t;
    edge_list_t edge_list_;

    node_id next_node_id_;

    Graph():
        next_node_id_(0)
    {}

    ~Graph()
    {
        // delete all nodes
        for (auto const & mapval : node_map_){
            node_ptr node = mapval.second;
            delete node;
        }
        // delete all edges
        for (auto edge : edge_list_){
            delete edge;
        }
    }

    size_t number_of_nodes() const { return node_map_.size(); }
    size_t number_of_edges() const { return edge_list_.size(); }

    /**
     * create a new node
     */
    node_ptr add_node(){
        node_ptr node = new Node(next_node_id_);
        next_node_id_++;
        node_map_.insert(std::pair<node_id, node_ptr> (node->id(), node));
        return node;
    }

    node_ptr add_node(node_id nodeid){
        node_ptr node = NULL;
        try {
            node = node_map_.at(nodeid);
            return node;
        } catch (std::out_of_range & e) {
            node = new Node(nodeid);
        }

        if (next_node_id_ < nodeid){
            next_node_id_ = nodeid + 1;
        }
        node_map_[node->id()] = node;
        return node;
    }


    /**
     * create a n new nodes
     */
    void add_nodes(node_id n){
        for (node_id i = 0; i < n; ++i){
            add_node();
        }
    }

    /**
     * return a pointer to the node with given node id
     */
    node_ptr get_node(node_id nodeid)
    {
        typedef std::map<node_id, node_ptr> maptype;
        maptype::iterator iter = node_map_.find(nodeid);
        if (iter == node_map_.end()){
            return NULL;
        }
        return iter->second;
        return NULL;
    }

    /**
     * add an edge from tail to head
     */
    edge_ptr add_edge(node_id tail, node_id head)
    {
        return _add_edge(get_node(tail), get_node(head));
    }
    edge_ptr _add_edge(node_ptr node_tail, node_ptr node_head)
    {
        assert(node_tail != NULL);
        assert(node_head != NULL);
        // check whether they're already connected
        edge_ptr edge = node_tail->get_successor_edge(node_head);
        if (edge != NULL){
            return edge;
        }
        edge = new Edge(node_tail, node_head);
        edge_list_.insert(edge);
        node_tail->add_out_edge(edge);
        node_head->add_in_edge(edge);
        return edge;
    }

    /**
     * remove a node and all edges connecting it
     */
    void remove_node(node_id nodeid){
        return _remove_node(get_node(nodeid));
    }

    void _remove_node(node_ptr u)
    {
        Node::edge_iterator eiter;

        // remove the edges from the nodes connected to u
        for (eiter = u->out_edge_begin(); eiter != u->out_edge_end(); ++eiter){
            edge_ptr edge = *eiter;
            edge->head()->remove_in_edge(edge);
        }
        for (eiter = u->in_edge_begin(); eiter != u->in_edge_end(); ++eiter){
            edge_ptr edge = *eiter;
            edge->tail()->remove_out_edge(edge);
        }

        Node::edge_list to_delete;
        to_delete.insert(u->in_edge_begin(), u->in_edge_end());
        to_delete.insert(u->out_edge_begin(), u->out_edge_end());

        // remove the edges from the edge list
        for (auto edge : to_delete){
            edge_list_.erase(edge);
        }

        // remove the node from the node list
        node_map_.erase(u->id());

        // deallocate the memory
        for (auto edge : to_delete){
            delete edge;
        }

        delete u;
    }

    /*
     * copy constructor
     */
    Graph(Graph & graph):
        next_node_id_(0)
    {
        for (auto const & mapval : graph.node_map_){
            node_ptr u = mapval.second;
//            std::cout << "making node " << u->id() << "\n";
            node_ptr unew = this->add_node(u->id());
            unew->tau = u->tau;
        }

        for (auto eiter = graph.edge_list_.begin(); eiter != graph.edge_list_.end(); ++eiter){
            edge_ptr edge = *eiter;
            node_id uid = edge->tail()->id();
            node_id vid = edge->head()->id();
            edge_ptr edge_new = this->add_edge(uid, vid);
            edge_new->P = edge->P;
        }
    }

};

inline std::ostream &operator<<(std::ostream &out, std::shared_ptr<Graph> g) {
    out << "nodes\n";
    out << "-----\n";
    for (auto const & nn : g->node_map_) {
        out << nn.first << " tau " << nn.second->tau << "\n";
    }
    out << "edges\n";
    out << "-----\n";
    for (auto const & e : g->edge_list_) {
        out << e->tail_->id() << " -> " << e->head_->id() << " P " << e->P << "\n";
    }
    return out;
}

}

#endif
