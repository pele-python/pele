#include <gtest/gtest.h>
#include "pele/graph.hpp"
#include "pele/ngt.hpp"

using pele::NGT;
using pele::node_id;

TEST(Graph, AddDuplicateEdges_NoEffect){
    pele::Graph g;
    g.add_node(0);
    g.add_node(1);
    ASSERT_EQ(g.number_of_edges(), 0u);
    g.add_edge(0, 1);
    ASSERT_EQ(g.number_of_edges(), 1u);
    g.add_edge(0, 1);
    ASSERT_EQ(g.number_of_edges(), 1u);
}

class NGT3 :  public ::testing::Test
{
public:
	NGT::rate_map_t rate_map;
	std::set<node_id> A, B;
    virtual void SetUp(){
        for (int i=0; i<3; ++i){
            for (int j=0; j<3; ++j){
                rate_map[std::pair<node_id, node_id>(i,j)] = 1;
            }
        }
        A.insert(0);
        B.insert(1);


    }
};

TEST_F(NGT3, Rates_Correct){
	NGT ngt(rate_map, A, B);
	ngt.compute_rates();
	ASSERT_NEAR(ngt.get_rate_AB(), 1, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA(), 1, 1e-9);
	ASSERT_NEAR(ngt.get_rate_AB_SS(), 1.5, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA_SS(), 1.5, 1e-9);
}

TEST_F(NGT3, RatesCommittors_Correct){
	NGT ngt(rate_map, A, B);
	ngt.compute_rates_and_committors();
	ASSERT_NEAR(ngt.get_rate_AB(), 1, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA(), 1, 1e-9);
	ASSERT_NEAR(ngt.get_rate_AB_SS(), 1.5, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA_SS(), 1.5, 1e-9);

	auto committors = ngt.get_committors();
	ASSERT_NEAR(committors[2], 0.5, 1e-9);
}



class NGT10 :  public ::testing::Test
{
public:
	NGT::rate_map_t rate_map;
	std::set<node_id> A, B;
	double kAB, kBA, kABSS, kBASS;
    virtual void SetUp(){
        for (int i=0; i<10; ++i){
            for (int j=0; j<10; ++j){
				if (i==j) continue;
                rate_map[std::pair<node_id, node_id>(i,j)] = double(i+j)/(i+1);
            }
        }
        A.insert(0);
        A.insert(1);
        A.insert(2);

        B.insert(3);
        B.insert(4);
        B.insert(5);
        kAB = 5.1013138820442565;
        kBA = 3;
        kABSS = 19.933950145409426;
        kBASS = 6.970856553435547;


    }
};

TEST_F(NGT10, Rates_Correct){
	NGT ngt(rate_map, A, B);
	ngt.compute_rates();
	ASSERT_NEAR(ngt.get_rate_AB(), kAB, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA(), kBA, 1e-9);
	ASSERT_NEAR(ngt.get_rate_AB_SS(), kABSS, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA_SS(), kBASS, 1e-9);
}

TEST_F(NGT10, RatesCommittors_Correct){
	NGT ngt(rate_map, A, B);
	ngt.compute_rates_and_committors();
	ASSERT_NEAR(ngt.get_rate_AB(), kAB, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA(), kBA, 1e-9);
	ASSERT_NEAR(ngt.get_rate_AB_SS(), kABSS, 1e-9);
	ASSERT_NEAR(ngt.get_rate_BA_SS(), kBASS, 1e-9);
}

