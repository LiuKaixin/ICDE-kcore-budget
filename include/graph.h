//
// Created by 刘凯鑫 on 2019-05-10.
//

#ifndef KCORELKX_GRAPH_H
#define KCORELKX_GRAPH_H

#include "mylib.h"
#include "config.h"

class Graph {
public:
    int n;
    long long m, kcore_size, kmb_core_size,km1_core_size;
    string data_folder;

    vector<int> k_tag,kc_nei;
    vector<vector<int>> adj_list,adj_list_copy;

    vector<int> core_num;
    vector<int> k_tag_tmp,nodes_less_k,del_nei;

    void init_nm();

    void convert_to_undirected_graph(string graph_path);

    Graph() = default;

    Graph(const string graph_path);

    void generate_subgraph_rand();

    void generate_subgraph_bfs();
    //void compute_kmb_kcore();
    void detect_kcore();
    void core_decompsition();
    int detect_kcore_with_anc(vector<int>& b);
    void save_subgraph(string outputFile);
};






#endif //KCORELKX_GRAPH_H
